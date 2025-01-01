import streamlit as st
import requests

BASE_API_URL = "http://127.0.0.1:8000/parameters"

def app():
    st.title("Parameters Management")

    st.subheader("Upload SLHA File")
    uploaded_file = st.file_uploader("Choose a SLHA file to upload", type=["slha", "lha"])
    if uploaded_file:
        response = requests.post(
            f"{BASE_API_URL}/upload",
            files={"file": (uploaded_file.name, uploaded_file.getvalue())}
        )
        if response.status_code == 200:
            st.success("File uploaded successfully!")
        else:
            st.error("Failed to upload file")

    st.subheader("Retrieve Parameter Value")
    block = st.text_input("Block Name")
    code = st.number_input("Parameter Code", step=1, min_value=0)
    if st.button("Get Parameter Value"):
        response = requests.get(f"{BASE_API_URL}/value", params={"block": block, "code": code})
        if response.status_code == 200:
            value = response.json().get("value")
            st.write(f"Value for block `{block}` and code `{code}`: {value}")
        else:
            st.error(response.json().get("detail", "Failed to retrieve parameter value"))
