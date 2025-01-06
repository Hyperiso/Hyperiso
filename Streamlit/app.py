import streamlit as st
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
page.app()