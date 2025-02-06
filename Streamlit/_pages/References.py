import streamlit as st
from Streamlit.Utils.common_elements import add_header, add_footer, apply_custom_background, apply_sidebar_style, apply_file_management_style
from Streamlit.Utils.common_elements import apply_custom_css

if "wide_mode" not in st.session_state:
    st.set_page_config(layout="wide", page_title="Hyperiso", page_icon="📊")
    st.session_state["wide_mode"] = True
    
def app():

    apply_custom_background()
    # add_header()
    apply_file_management_style()
    apply_custom_css()
    apply_sidebar_style()
    st.title("References")
    st.write("1. FastAPI Documentation: [FastAPI](https://fastapi.tiangolo.com)")
    st.write("2. Streamlit Documentation: [Streamlit](https://docs.streamlit.io)")
    st.write("3. Passlib Documentation: [Passlib](https://passlib.readthedocs.io)")
    st.write("4. JOSE Documentation: [PyJWT](https://pyjwt.readthedocs.io)")

    add_footer()

app()