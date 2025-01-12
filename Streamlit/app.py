import streamlit as st
st.set_page_config(layout="wide", page_title="Hyperiso", page_icon="📊")
from pages import login, parameters, wilson, observables

PAGES = {
    "Login": login,
    "Parameters": parameters,
    "Wilson Coefficients": wilson,
    "Observables": observables,
}

st.sidebar.title("Navigation")
selection = st.sidebar.radio("Go to", list(PAGES.keys()))
page = PAGES[selection]
if hasattr(page, "app"):
    page.app()
else:
    st.error(f"La page sélectionnée n'a pas de fonction `app` : {selection}")