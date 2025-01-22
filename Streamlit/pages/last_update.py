import streamlit as st
from Streamlit.Utils.common_elements import add_header, add_footer, apply_custom_background, apply_sidebar_style, apply_file_management_style
from Streamlit.Utils.common_elements import apply_custom_css

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

if __name__ == "__main__":
    app()