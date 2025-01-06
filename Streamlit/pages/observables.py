import streamlit as st
import requests

BASE_API_URL = "http://127.0.0.1:8000"

def app():
    st.title("Observable Calculations")

    st.subheader("Calculate an Observable")
    observable_name = st.text_input("Observable Name")
    if st.button("Calculate Observable"):
        response = requests.get(f"{BASE_API_URL}/observables/calculate", params={"name": observable_name})
        if response.status_code == 200:
            result = response.json().get("value")
            st.write(f"Value of `{observable_name}`: {result}")
        else:
            st.error("Failed to calculate observable")

    st.subheader("Plot Observable Variation")
    param_name = st.text_input("Parameter Name for Variation")
    param_min = st.number_input("Minimum Value", step=0.1)
    param_max = st.number_input("Maximum Value", step=0.1)
    if st.button("Plot Observable Variation"):
        response = requests.get(
            f"{BASE_API_URL}/observables/plot_variation",
            params={"param_name": param_name, "min_value": param_min, "max_value": param_max}
        )
        if response.status_code == 200:
            plot_data = response.json()
            st.line_chart(plot_data["values"])
        else:
            st.error("Failed to generate plot")

    st.subheader("Calculate Chi2")
    if st.button("Calculate Chi2"):
        response = requests.get(f"{BASE_API_URL}/observables/chi2")
        if response.status_code == 200:
            chi2_value = response.json().get("value")
            st.write(f"Chi2 Value: {chi2_value}")
        else:
            st.error("Failed to calculate Chi2")
