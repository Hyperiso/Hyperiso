import streamlit as st
import requests
from Streamlit.Utils.common_elements import add_header, add_footer, apply_custom_background, apply_sidebar_style, apply_file_management_style
from Streamlit.Utils.common_elements import apply_custom_css
import os
import matplotlib.pyplot as plt
import numpy as np

BASE_API_URL = "http://127.0.0.1:8000"

if "selected_file" not in st.session_state:
    st.session_state.selected_file = None

if "selected_model" not in st.session_state:
    st.session_state.selected_model = "SM"

if "param_type" not in st.session_state:
    st.session_state.param_type = "SM"

if "use_marty" not in st.session_state:
    st.session_state.use_marty = False

if "is_spectrum" not in st.session_state:
    st.session_state.is_spectrum = False

if "has_wilson" not in st.session_state:
    st.session_state.has_wilson = False

if "has_obs" not in st.session_state:
    st.session_state.has_obs = False

if "group" not in st.session_state:
    st.session_state.group = "BCoefficientGroup"

if "order" not in st.session_state:
    st.session_state.order = "LO"

if "coeff_name" not in st.session_state:
    st.session_state.coeff_name = "C1"

if "base" not in st.session_state:
    st.session_state.base = "base1"

if "matching_scale" not in st.session_state:
    st.session_state.matching_scale = 81

groups = ["BCoefficientGroup", "BPrimeCoefficientGroup", "BScalarCoefficientGroup"]

coeff_by_group = {
    "BCoefficientGroup" : ["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"],
    "BPrimeCoefficientGroup" : ["CP1", "CP2", "CP3", "CP4", "CP5", "CP6", "CP7", "CP8", "CP9", "CP10", "CPQ1", "CPQ2"],
    "BScalarCoefficientGroup" : ["CQ1", "CQ2"]
}

running_base_b = ["base1", "base2"]

def truc():
    print(st.session_state.use_marty, st.session_state.selected_file, st.session_state.selected_model)

def app():
    apply_custom_background()
    add_header()
    apply_file_management_style()
    apply_custom_css()
    st.title("Wilson Coefficients")

    st.sidebar.header("Configuration")
    model = st.sidebar.selectbox("Select Model", ["SM", "SUSY", "THDM", "CUSTOM"], key="selected_model")
    lha_files = os.listdir("DataBase/lha")
    # lha_files = ["file1.lha", "file2.lha", "file3.lha"]  # Simuler une liste de fichiers LHA
    lha_file = st.sidebar.selectbox("Select LHA File", lha_files, key="selected_file")
    use_marty = st.sidebar.checkbox("Use Marty (forces LO Order)", key="use_marty")

    if st.sidebar.button("Set LHA and Model"):
        if st.session_state.selected_file != None:
            response = requests.post(
                f"{BASE_API_URL}/parameters/set_memory_manager",
                json={"lha_file": os.path.join("DataBase/lha",st.session_state.selected_file), "model": st.session_state.selected_model,
                    "use_marty": st.session_state.use_marty,
                "is_spectrum": st.session_state.is_spectrum,
                "has_wilsons": st.session_state.has_wilson,
                "has_obs": st.session_state.has_obs}
            )
            if response.status_code == 200:
                st.sidebar.success("LHA and Model configured!")
            else:
                st.sidebar.error("Failed to configure LHA and Model")

    st.sidebar.subheader("Global Parameters")
    group = st.sidebar.selectbox("CoefficientGroup", groups, key = "group")
    matching_scale = st.sidebar.number_input("Q_match (Matching Scale)", min_value=81.0, step=0.1, key = "matching_scale")
    q_value = st.sidebar.number_input("Q (Running Scale)", min_value=1.0, step=0.1)
    qcd_order = st.sidebar.selectbox("QCD Order", ["LO", "NLO", "NNLO"], disabled=use_marty, key = "qcd_order")
    base = st.sidebar.selectbox("running base" , running_base_b, key = "base")

    if st.sidebar.button("Set Global Parameters"):
        register_group = requests.post(f"{BASE_API_URL}/wilson/register_group", json = {"group" : group, "model" : st.session_state.selected_model})
        set_q_match = requests.post(f"{BASE_API_URL}/wilson/set_matching_scale", json={"model" : st.session_state.selected_model, "group" : group, "scale": matching_scale, "qcd_order" : qcd_order})
        set_group_scale = requests.post(f"{BASE_API_URL}/wilson/set_group_scale", json={"model" : st.session_state.selected_model, "group" : group, "scale": q_value, "qcd_order" : qcd_order})

        if set_q_match.status_code == 200 and set_group_scale.status_code == 200:
            st.sidebar.success("Global parameters configured!")
        else:
            st.sidebar.error("Failed to configure global parameters")

    col1, col2, col3 = st.columns(3)

    with col1:
        st.subheader("Matching Analysis")
        # group = st.selectbox("Select Group", ["Group1", "Group2", "Group3"])
        coeff_name = st.selectbox("Select Coefficient", coeff_by_group[st.session_state.group]+["all"], key="coeff_name")
        # coeff_name = st.text_input("Coefficient Name")
        # coeff_order = st.selectbox("Order", ["LO", "NLO", "NNLO"])
        if st.button("Retrieve Matching Coefficient"):
            if coeff_name == "all":
                st.write(f"Group: {group}")
                for coef in coeff_by_group[st.session_state.group]:
                    response = requests.get(
                        f"{BASE_API_URL}/wilson/get_coefficient",
                        params={"model" : st.session_state.selected_model, "group": st.session_state.group, "name": coef, "order": st.session_state.order},
                    )
                    if response.status_code == 200:
                        data = response.json()
                        st.write(f"Coefficient: {coef}")
                        st.write(f"Value: {round(data['coeff_real'],5)} + i{round(data['coeff_img'], 5)}")
                    else:
                        st.error("Failed to retrieve coefficient")

            elif coeff_name != "all":
                response = requests.get(
                    f"{BASE_API_URL}/wilson/get_coefficient",
                    params={"model" : st.session_state.selected_model, "group": st.session_state.group, "name": st.session_state.coeff_name, "order": st.session_state.order},
                )
                if response.status_code == 200:
                    data = response.json()
                    st.write(f"Group: {group}")
                    st.write(f"Coefficient: {coeff_name}")
                    st.write(f"Value: {round(data['coeff_real'],5)} + i{round(data['coeff_img'], 5)}")
                else:
                    st.error("Failed to retrieve coefficient")

        st.subheader("Coefficient Variation")
        param_block = st.text_input("Parameter Block")
        param_pdgcode = st.number_input("Parameter pdgCode", step =1)
        param_min = st.number_input("Min Value", step=0.1)
        param_max = st.number_input("Max Value", step=0.1)
        param_steps = st.number_input("Number of point", min_value=2, step=1, value=10)
        if st.button("Plot Parameter Variation"):
            response = requests.get(
                f"{BASE_API_URL}/wilson/plot_coefficients",
                json={
                    "group": st.session_state.group,
                    "name": st.session_state.coeff_name,
                    "order": st.session_state.order,
                    "matching_scale" : st.session_state.matching_scale,
                    "param_block": param_block,
                    "param_pdgcode" : param_pdgcode,
                    "min_value": param_min,
                    "max_value": param_max,
                    "steps": param_steps,
                },
            )
            if response.status_code == 200:
                plot_data = response.json()
                fig, ax = plt.subplots()
                plt.plot(np.linspace(param_min, param_max, param_steps), [val["coefficient"] for val in plot_data["values"]])
                st.pyplot(fig)
            else:
                st.error("Failed to generate plot")

    with col2:
        st.subheader("Running Analysis")
        if st.button("Retrieve Coefficients by Q"):
            response = requests.get(
                f"{BASE_API_URL}/calculate_coefficient",
                params={"group": group, "name": coeff_name, "order": coeff_order},
            )
            if response.status_code == 200:
                st.write(f"Value at Q={q_value}: {response.json()['value']}")
            else:
                st.error("Failed to retrieve coefficient")

    with col3:
        st.subheader("Analyze All LHAs")
        if st.button("Compute for All LHAs"):
            response = requests.get(
                f"{BASE_API_URL}/calculate_all_lhas",
                params={"group": group, "name": coeff_name, "order": coeff_order},
            )
            if response.status_code == 200:
                lha_data = response.json()["results"]
                st.bar_chart([result["value"] for result in lha_data])
            else:
                st.error("Failed to compute for all LHAs")
    add_footer()

if __name__ == "__main__":
    app()