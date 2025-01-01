import streamlit as st
import requests

BASE_API_URL = "http://127.0.0.1:8000"

def app():
    st.title("Wilson Coefficients")

    st.subheader("Set Matching Scale")
    matching_scale = st.number_input("Matching Scale (GeV)", min_value=0.0, step=0.1)
    if st.button("Set Matching Scale"):
        response = requests.post(
            f"{BASE_API_URL}/wilson/set_matching_scale", json={"scale": matching_scale}
        )
        if response.status_code == 200:
            st.success("Matching scale set successfully!")
        else:
            st.error("Failed to set matching scale")

    st.subheader("Retrieve Wilson Coefficient")
    coefficient_name = st.text_input("Coefficient Name")
    if st.button("Get Coefficient Values"):
        response = requests.get(f"{BASE_API_URL}/wilson/coefficient", params={"name": coefficient_name})
        if response.status_code == 200:
            data = response.json()
            st.write(f"Matching Value: {data['matching_value']}")
            st.write(f"Running Value: {data['running_value']}")
        else:
            st.error("Failed to retrieve coefficient values")

    st.subheader("Plot Coefficient Variation")
    param_name = st.text_input("Parameter Name for Variation")
    param_min = st.number_input("Minimum Value", step=0.1)
    param_max = st.number_input("Maximum Value", step=0.1)
    if st.button("Plot Coefficient Variation"):
        response = requests.get(
            f"{BASE_API_URL}/wilson/plot_variation",
            params={"param_name": param_name, "min_value": param_min, "max_value": param_max}
        )
        if response.status_code == 200:
            plot_data = response.json()
            st.line_chart(plot_data["values"])
        else:
            st.error("Failed to generate plot")
