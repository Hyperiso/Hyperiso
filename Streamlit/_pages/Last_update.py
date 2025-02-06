import streamlit as st
from Streamlit.Utils.common_elements import add_header, add_footer, apply_custom_background, apply_sidebar_style, apply_file_management_style
from Streamlit.Utils.common_elements import apply_custom_css

if "wide_mode" not in st.session_state:
    st.set_page_config(layout="wide", page_title="Hyperiso", page_icon="📊")
    st.session_state["wide_mode"] = True
    
def app():

    apply_custom_background()
    add_header()
    apply_file_management_style()
    apply_custom_css()
    apply_sidebar_style()
    st.title("Last Update")
    st.write("Dernière mise à jour : 20 janvier 2025")
    st.write("Cette application a été mise à jour pour inclure des fonctionnalités de sécurité supplémentaires.")
    add_footer()


app()