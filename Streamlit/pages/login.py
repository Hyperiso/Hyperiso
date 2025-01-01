import streamlit as st
import requests

def app():
    st.title("Login")

    username = st.text_input("Username")
    password = st.text_input("Password", type="password")
    if st.button("Login"):
        response = requests.post(
            "http://127.0.0.1:8000/auth/login",
            data={"username": username, "password": password},
        )
        if response.status_code == 200:
            st.session_state["token"] = response.json()["access_token"]
            st.success("Logged in successfully")
        else:
            st.error("Invalid credentials")
