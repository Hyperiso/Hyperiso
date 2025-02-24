import streamlit as st
import requests

from Streamlit.Utils.common_elements import add_header, add_footer, apply_custom_background, apply_sidebar_style, apply_file_management_style
from Streamlit.Utils.common_elements import apply_custom_css

if "wide_mode" not in st.session_state:
    st.set_page_config(layout="wide", page_title="Hyperiso", page_icon="📊")
    st.session_state["wide_mode"] = True
    
BASE_API_URL = "http://127.0.0.1:8000"

if "selected_model" not in st.session_state:
    st.session_state.selected_model = "SM"


if "obs_list" not in st.session_state:
    _ = requests.get(f"{BASE_API_URL}/observables/get_obs_names").json()
    st.session_state.obs_list = _["obs_names"]

def app():
    apply_custom_background()
    apply_file_management_style()
    apply_custom_css()
    apply_sidebar_style()
    st.title(r"Observable and $\chi^2$ Calculations")

    obs_type = st.selectbox("What kind of analysis you want ?", ["Single Observable Analysis", "Uncertainties calculation", r"$\Chi^2 calculation"], key = "obs_type")

    if obs_type == "Single Observable Analysis":
        obs_name = st.selectbox("Select Observable", st.session_state.obs_list, key="obs_name")
        obs_order = st.selectbox("Select QCDOrder", ["LO", "NLO", "NNLO"])
        if st.button("Calculate Observable"):
            print(obs_name, obs_order)
            _ = requests.post(f"{BASE_API_URL}/observables/add_observable", json={"obs": obs_name, "order" : obs_order})
            print(_.json())
            response = requests.get(f"{BASE_API_URL}/observables/compute_observable", params={"name": obs_name})
            if response.status_code == 200:
                result = response.json().get("value")
                st.write(f"Value of `{obs_name}`: {result}")
            else:
                st.error("Failed to calculate observable")

    if obs_type == "Uncertainties calculation":
        obs_names = st.multiselect("Select Observables", st.session_state.obs_list, key="obs_names")
        b = 1

    if obs_type == r"$\Chi^2 calculation":
        obs_names_chi = st.checkbox("Select Observables", st.session_state.obs_list, key="obs_names_chi")
        c = 1

    # st.subheader("Calculate an Observable")
    # observable_name = st.text_input("Observable Name")
    # if st.button("Calculate Observable"):
    #     response = requests.get(f"{BASE_API_URL}/observables/compute_observable", params={"name": observable_name})
    #     if response.status_code == 200:
    #         result = response.json().get("value")
    #         st.write(f"Value of `{observable_name}`: {result}")
    #     else:
    #         st.error("Failed to calculate observable")

    # st.subheader("Plot Observable Variation")
    # param_name = st.text_input("Parameter Name for Variation")
    # param_min = st.number_input("Minimum Value", step=0.1)
    # param_max = st.number_input("Maximum Value", step=0.1)
    # if st.button("Plot Observable Variation"):
    #     response = requests.get(
    #         f"{BASE_API_URL}/observables/plot_variation",
    #         params={"param_name": param_name, "min_value": param_min, "max_value": param_max}
    #     )
    #     if response.status_code == 200:
    #         plot_data = response.json()
    #         st.line_chart(plot_data["values"])
    #     else:
    #         st.error("Failed to generate plot")

    # st.subheader("Calculate Chi2")
    # if st.button("Calculate Chi2"):
    #     response = requests.get(f"{BASE_API_URL}/observables/chi2")
    #     if response.status_code == 200:
    #         chi2_value = response.json().get("value")
    #         st.write(f"Chi2 Value: {chi2_value}")
    #     else:
    #         st.error("Failed to calculate Chi2")
    add_footer()

app()