import streamlit as st
st.set_page_config(layout="wide", page_title="Hyperiso", page_icon="📊")
from pages import login, generation, parameters, wilson, observables, references, last_update

PAGES = {
    "Login": login,
    "LHA Generation" : generation,
    "Parameters": parameters,
    "Wilson Coefficients": wilson,
    "Observables": observables,
    "Last update": last_update,
    "References": references
}

if "authenticated" not in st.session_state:
    st.session_state["authenticated"] = False

st.sidebar.title("Navigation")
if not st.session_state["authenticated"]:
    st.warning("Veuillez vous connecter pour accéder à l'application.")
    login.app()
else:
    selection = st.sidebar.radio("Go to", list(PAGES.keys()))
    page = PAGES[selection]
    if hasattr(page, "app"):
        page.app()
    else:
        st.error(f"La page sélectionnée n'a pas de fonction `app` : {selection}")
