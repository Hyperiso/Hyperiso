import streamlit as st
import requests
from Streamlit.Utils.common_elements import add_header, add_footer, apply_custom_background, apply_sidebar_style, apply_file_management_style
from Streamlit.Utils.common_elements import apply_custom_css

def app():
    st.title("Login")
    apply_custom_background()
    add_header()
    apply_file_management_style()
    apply_custom_css()
    apply_sidebar_style()
    
    if "authenticated" not in st.session_state:
        st.session_state["authenticated"] = False
    
    if st.session_state["authenticated"]:
        st.success("Vous êtes déjà connecté.")
        st.stop()
    
    username = st.text_input("Username")
    password = st.text_input("Password", type="password")
    if st.button("Login"):
        response = requests.post(
            "http://127.0.0.1:8000/auth/login",
            data={"username": username, "password": password},
        )
        if response.status_code == 200:
            st.session_state["authenticated"] = True
            st.session_state["token"] = response.json()["access_token"]
            st.success("Logged in successfully")
        else:
            st.error("Invalid credentials")
    
    add_footer()

if __name__ == "__main__":
    app()