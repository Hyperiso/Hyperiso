# dash_bot/data/sources.py
import pandas as pd
from typing import List
from DataAnalysis.config import DATA_PATH

try:
    import yfinance as yf  # type: ignore
    HAS_YFINANCE = True
except Exception:
    HAS_YFINANCE = False

try:
    import ccxt  # type: ignore
    HAS_CCXT = True
except Exception:
    HAS_CCXT = False

def _map_usdt_to_yf(symbol_usdt: str) -> str:
    # 'BTCUSDT' -> 'BTC-USD'
    s = symbol_usdt.upper()
    if s.endswith("USDT"):
        base = s[:-4]
        return f"{base}-USD"
    return s

def fetch_prices_yfinance(symbols: List[str], start: str, end: str, interval: str) -> pd.DataFrame:
    if not HAS_YFINANCE:
        raise RuntimeError("yfinance non installé: pip install yfinance")
    tickers = [_map_usdt_to_yf(s) for s in symbols]
    df = yf.download(
        tickers=tickers,
        start=start,
        end=end,
        interval=interval,
        auto_adjust=False,
        progress=False,
        group_by="ticker",
        threads=True,
    )
    out = {}
    if isinstance(df.columns, pd.MultiIndex):
        for t, s in zip(tickers, symbols):
            if (t, "Close") in df.columns:
                out[s] = df[(t, "Close")]
    else:
        # single ticker
        if "Close" in df.columns:
            out[symbols[0]] = df["Close"]

    res = pd.DataFrame(out).sort_index()
    res.index = pd.to_datetime(res.index)
    return res

def fetch_prices_binance_ccxt(symbols, *, timeframe: str, lookback_days: int, limit: int | None = None):
    if not HAS_CCXT:
        raise RuntimeError("ccxt non installé: pip install ccxt")
    # Binance futures: binanceusdm (USDT-margined perpetual futures)
    ex = ccxt.binanceusdm({"enableRateLimit": True})
    ex.load_markets()

    tf_map = {"5m": "5m", "15m": "15m", "1h": "1h", "4h": "4h", "1d": "1d"}
    timeframe = tf_map.get(interval, interval)

    since_ms = int(pd.to_datetime(start).timestamp() * 1000)
    end_ms = int(pd.to_datetime(end).timestamp() * 1000)

    out: Dict[str, pd.Series] = {}
    errors: Dict[str, str] = {}

    for sym in symbols:
        market_sym = sym.replace("USDT", "/USDT")
        if market_sym not in ex.markets:
            errors[sym] = f"Market not found on binanceusdm: {market_sym}"
            continue

        all_rows = []
        since = since_ms
        try:
            # paginate
            while True:
                ohlcv = ex.fetch_ohlcv(market_sym, timeframe=timeframe, since=since, limit=1500)
                if not ohlcv:
                    break
                all_rows.extend(ohlcv)
                last = ohlcv[-1][0]
                since = last + 1
                if last >= end_ms or len(ohlcv) < 1500:
                    break
        except Exception as e:
            errors[sym] = f"fetch_ohlcv error: {type(e).__name__}: {e}"
            continue

        if not all_rows:
            errors[sym] = "No OHLCV returned (empty)."
            continue

        arr = np.array(all_rows, dtype=float)
        idx = pd.to_datetime(arr[:, 0], unit="ms")
        close = pd.Series(arr[:, 4], index=idx)
        close = close[(close.index >= pd.to_datetime(start)) & (close.index < pd.to_datetime(end))]
        if close.dropna().empty:
            errors[sym] = "Only NaNs after date filter."
            continue
        out[sym] = close

    df = pd.DataFrame(out).sort_index()
    return df, errors

def load_prices_csv(data_path: str, *, tz: str | None = None):
    from DataAnalysis.src.data_fetching import fetch_price_data
    from DataAnalysis.Utils.load_prices import load_all_prices
    fetch_price_data()
    df = load_all_prices(data_path)
    df.index = pd.to_datetime(df.index)
    return df.sort_index()
