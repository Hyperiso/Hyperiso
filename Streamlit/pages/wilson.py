import streamlit as st
import requests
from Streamlit.Utils.common_elements import add_header, add_footer, apply_custom_background, apply_sidebar_style, apply_file_management_style
from Streamlit.Utils.common_elements import apply_custom_css
import os
import matplotlib.pyplot as plt
import numpy as np

BASE_API_URL = "http://127.0.0.1:8000"

def coef_name_transformation(coeff_name : str) -> str:
    for i in range(coeff_name.__len__()):
        if coeff_name[i:].isdigit():
            return rf"${coeff_name[0:i]}_{coeff_name[i:]}$"

# Initialisation des états de session
if "selected_file" not in st.session_state:
    st.session_state.selected_file = None

if "selected_model" not in st.session_state:
    st.session_state.selected_model = "SM"

if "is_configured" not in st.session_state:
    st.session_state.is_configured = False

if "global_params_set" not in st.session_state:
    st.session_state.global_params_set = False

groups = ["BCoefficientGroup", "BPrimeCoefficientGroup", "BScalarCoefficientGroup"]

coeff_by_group = {
    "BCoefficientGroup": ["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"],
    "BPrimeCoefficientGroup": ["CP1", "CP2", "CP3", "CP4", "CP5", "CP6", "CP7", "CP8", "CP9", "CP10", "CPQ1", "CPQ2"],
    "BScalarCoefficientGroup": ["CQ1", "CQ2"]
}

running_base_b = ["base1", "base2"]

def app():
    apply_custom_background()
    apply_sidebar_style()
    add_header()
    apply_file_management_style()
    apply_custom_css()
    st.title("Wilson Coefficients")

    st.sidebar.header("Configuration")
    model = st.sidebar.selectbox("Select Model", ["SM", "SUSY", "THDM", "CUSTOM"], key="selected_model")
    lha_files = os.listdir("DataBase/lha")
    lha_file = st.sidebar.selectbox("Select LHA File", lha_files, key="selected_file")
    use_marty = st.sidebar.checkbox("Use Marty (forces LO Order)", key="use_marty")

    if st.sidebar.button("Set LHA and Model"):
        if st.session_state.selected_file:
            response = requests.post(
                f"{BASE_API_URL}/parameters/set_memory_manager",
                json={
                    "lha_file": os.path.join("DataBase/lha", st.session_state.selected_file),
                    "model": st.session_state.selected_model,
                    "use_marty": st.session_state.use_marty,
                },
            )
            if response.status_code == 200:
                st.sidebar.success("LHA and Model configured!")
                st.session_state.is_configured = True
            else:
                st.sidebar.error("Failed to configure LHA and Model")

    if st.session_state.is_configured:
        st.sidebar.subheader("Global Parameters")
        group = st.sidebar.selectbox("CoefficientGroup", groups, key="group", disabled=not st.session_state.is_configured)
        matching_scale = st.sidebar.number_input("Q_match (Matching Scale)", min_value=81.0, step=0.1, key="matching_scale")
        running_scale = st.sidebar.number_input("Q (Running Scale)", min_value=1.0, step=0.1, key = "running_scale")
        qcd_order = st.sidebar.selectbox("QCD Order", ["LO", "NLO", "NNLO"], key="qcd_order")
        base = st.sidebar.selectbox("Running base", running_base_b, key="base")

        if st.sidebar.button("Set Global Parameters"):
            register_group = requests.post(f"{BASE_API_URL}/wilson/register_group", json={"group": group, "model": model})
            set_q_match = requests.post(
                f"{BASE_API_URL}/wilson/set_matching_scale",
                json={"model": model, "group": group, "scale": matching_scale, "qcd_order": qcd_order},
            )
            set_group_scale = requests.post(
                f"{BASE_API_URL}/wilson/set_group_scale",
                json={"model": model, "group": group, "scale": running_scale, "qcd_order": qcd_order},
            )

            if set_q_match.status_code == 200 and set_group_scale.status_code == 200:
                st.sidebar.success("Global parameters configured!")
                st.session_state.global_params_set = True
            else:
                st.sidebar.error("Failed to configure global parameters")

    if st.session_state.global_params_set:
        st.subheader("Coefficient Selection")
        coeff_name = st.selectbox("Select Coefficient", coeff_by_group[st.session_state.group] + ["all"], key="coeff_name")

        col1, col2, col3 = st.columns(3)

        with col1:
            st.subheader("Matching Analysis")
            # group = st.selectbox("Select Group", ["Group1", "Group2", "Group3"])
            # coeff_name = st.selectbox("Select Coefficient", coeff_by_group[st.session_state.group]+["all"], key="coeff_name")
            # coeff_name = st.text_input("Coefficient Name")
            # coeff_order = st.selectbox("Order", ["LO", "NLO", "NNLO"])
            if st.button("Retrieve Matching Coefficient"):
                if coeff_name == "all":
                    st.write(f"Group: {group}")
                    for coef in coeff_by_group[st.session_state.group]:
                        response = requests.get(
                            f"{BASE_API_URL}/wilson/get_coefficient",
                            params={"model" : st.session_state.selected_model, "group": st.session_state.group, "name": coef, "order": st.session_state.qcd_order},
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
                        json={"model" : st.session_state.selected_model, "group": st.session_state.group, "name": st.session_state.coeff_name, "qcd_order": st.session_state.qcd_order, "scale" : st.session_state.matching_scale},
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
                        "order": st.session_state.qcd_order,
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
                    plt.grid(True)
                    plt.tick_params("both", which = "major", direction="in", length = 5)
                    plt.tick_params("both", which = "minor", direction="in", length = 2.5)
                    plt.xlabel(f"{param_block}({param_pdgcode})")
                    plt.ylabel(coef_name_transformation(st.session_state.coeff_name))
                    plt.xlim(param_min, param_max)
                    st.pyplot(fig)
                else:
                    st.error("Failed to generate plot")
                
                response = requests.post(
                    f"{BASE_API_URL}/parameters/set_memory_manager",
                    json={
                        "lha_file": os.path.join("DataBase/lha", st.session_state.selected_file),
                        "model": st.session_state.selected_model,
                        "use_marty": st.session_state.use_marty,
                    },
                )

        with col2:
            st.subheader("Running Analysis")
            # group = st.selectbox("Select Group", ["Group1", "Group2", "Group3"])
            # coeff_name = st.selectbox("Select Coefficient", coeff_by_group[st.session_state.group]+["all"], key="coeff_name")
            # coeff_name = st.text_input("Coefficient Name")
            # coeff_order = st.selectbox("Order", ["LO", "NLO", "NNLO"])
            if st.button("Retrieve running Coefficient"):
                if coeff_name == "all":
                    st.write(f"Group: {group}")
                    for coef in coeff_by_group[st.session_state.group]:
                        response = requests.get(
                            f"{BASE_API_URL}/wilson/get_run_coefficient",
                            params={"model" : st.session_state.selected_model, "group": st.session_state.group, "name": coef, "order": st.session_state.qcd_order},
                        )
                        if response.status_code == 200:
                            data = response.json()
                            st.write(f"Coefficient: {coef}")
                            st.write(f"Value: {round(data['coeff_real'],5)} + i{round(data['coeff_img'], 5)}")
                        else:
                            st.error("Failed to retrieve coefficient")

                elif coeff_name != "all":
                    response = requests.get(
                        f"{BASE_API_URL}/wilson/get_run_coefficient",
                        params={"model" : st.session_state.selected_model, "group": st.session_state.group, "name": st.session_state.coeff_name, "order": st.session_state.qcd_order},
                    )
                    if response.status_code == 200:
                        data = response.json()
                        st.write(f"Group: {group}")
                        st.write(f"Coefficient: {coeff_name}")
                        st.write(f"Value: {round(data['coeff_real'],5)} + i{round(data['coeff_img'], 5)}")
                    else:
                        st.error("Failed to retrieve coefficient")

            st.subheader("Coefficient Variation given a Parameters")
            param_block = st.text_input("Parameter Block", key="param_block_run")
            param_pdgcode = st.number_input("Parameter pdgCode", step =1, key="param_pdgcode_run")
            param_min = st.number_input("Min Value", step=0.1, key="min_value_run")
            param_max = st.number_input("Max Value", step=0.1, key="max_value_run")
            param_steps = st.number_input("Number of point", min_value = 2, step=1, value=10, key = "step_value_run")
    
            if st.button("Plot Running wilson coefficient Variation"):
                response = requests.get(
                    f"{BASE_API_URL}/wilson/plot_run_coefficients",
                    json={
                        "group": st.session_state.group,
                        "name": st.session_state.coeff_name,
                        "order": st.session_state.qcd_order,
                        "matching_scale" : st.session_state.matching_scale,
                        "running_scale" : st.session_state.running_scale,
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
                    plt.plot(np.linspace(param_min, param_max, param_steps), [val["coefficient"] for val in plot_data["values"]], color = "red")
                    plt.grid(True)
                    plt.tick_params("both", which = "major", direction="in", length = 5)
                    plt.tick_params("both", which = "minor", direction="in", length = 2.5)
                    plt.xlabel(f"{param_block}({param_pdgcode})")
                    plt.ylabel(coef_name_transformation(st.session_state.coeff_name))
                    plt.xlim(param_min, param_max)
                    st.pyplot(fig)
                else:
                    st.error("Failed to generate plot")
                response = requests.post(
                        f"{BASE_API_URL}/parameters/set_memory_manager",
                        json={
                            "lha_file": os.path.join("DataBase/lha", st.session_state.selected_file),
                            "model": st.session_state.selected_model,
                            "use_marty": st.session_state.use_marty,
                        },
                    )
            st.subheader("Coefficient Variation given running scale")

            Q_min = st.number_input("Min Q Value", step=0.1, key="min_Q_run")
            Q_max = st.number_input("Max Q Value", step=0.1, key="max_Q_run")
            Q_step = st.number_input("Number of point", min_value = 2, step=1, value=10, key = "step_Q_run")
    
            if st.button("Plot Running wilson coefficient Variation with Q"):
                response = requests.get(
                    f"{BASE_API_URL}/wilson/plot_run_coefficients_Q",
                    json={
                        "group": st.session_state.group,
                        "name": st.session_state.coeff_name,
                        "order": st.session_state.qcd_order,
                        "matching_scale" : st.session_state.matching_scale,
                        "min_value": Q_min,
                        "max_value": Q_max,
                        "steps": Q_step,
                    },
                )
                if response.status_code == 200:
                    plot_data = response.json()
                    fig, ax = plt.subplots()
                    plt.plot(np.linspace(Q_min, Q_max, Q_step), [val["coefficient"] for val in plot_data["values"]], color = "red")
                    plt.grid(True)
                    plt.tick_params("both", which = "major", direction="in", length = 5)
                    plt.tick_params("both", which = "minor", direction="in", length = 2.5)
                    plt.xlabel(f"Q [GeV]")
                    plt.ylabel(coef_name_transformation(st.session_state.coeff_name))
                    plt.xlim(Q_min, Q_max)
                    st.pyplot(fig)
                else:
                    st.error("Failed to generate plot")

        with col3:
            st.subheader("Analyze All LHAs")
            if st.button("Compute for All LHAs"):
                response = requests.get(
                    f"{BASE_API_URL}/wilson/calculate_all_lhas",
                    params={"group": group, "name": coeff_name, "order": st.session_state.qcd_order, "scale" : st.session_state.matching_scale},
                )
                if response.status_code == 200:
                    lha_data = response.json()["results"]
                    st.bar_chart([result["value"] for result in lha_data])
                else:
                    st.error("Failed to compute for all LHAs")
    add_footer()

if __name__ == "__main__":
    app()