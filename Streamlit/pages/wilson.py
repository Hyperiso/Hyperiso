import streamlit as st
import requests
from Streamlit.Utils.common_elements import add_header, add_footer, apply_custom_background, apply_sidebar_style, apply_file_management_style


BASE_API_URL = "http://127.0.0.1:8000"

def app():
    apply_custom_background()
    add_header() 
    apply_file_management_style()
    st.title("Wilson Coefficients")

    st.sidebar.header("Configuration")
    model = st.sidebar.selectbox("Select Model", ["SM", "SUSY", "THDM", "CUSTOM"])
    lha_files = ["file1.lha", "file2.lha", "file3.lha"]  # Simuler une liste de fichiers LHA
    lha_file = st.sidebar.selectbox("Select LHA File", lha_files)
    use_marty = st.sidebar.checkbox("Use Marty (forces LO Order)")

    if st.sidebar.button("Set LHA and Model"):
        response = requests.post(
            f"{BASE_API_URL}/set_lha_model",
            json={
                "lha_file": lha_file,
                "model": model,
                "use_marty": use_marty,
            },
        )
        if response.status_code == 200:
            st.sidebar.success("LHA and Model configured!")
        else:
            st.sidebar.error("Failed to configure LHA and Model")

    st.sidebar.subheader("Global Parameters")
    q_match = st.sidebar.number_input("Q_match (Matching Scale)", min_value=0.0, step=0.1)
    q_value = st.sidebar.number_input("Q (Running Scale)", min_value=0.0, step=0.1)
    qcd_order = st.sidebar.selectbox("QCD Order", ["LO", "NLO", "NNLO"], disabled=use_marty)

    if st.sidebar.button("Set Global Parameters"):
        set_q_match = requests.post(f"{BASE_API_URL}/wilson/set_matching_scale", json={"scale": q_match})
        set_group_scale = requests.post(f"{BASE_API_URL}/wilson/set_group_scale", json={"scale": q_value})

        if set_q_match.status_code == 200 and set_group_scale.status_code == 200:
            st.sidebar.success("Global parameters configured!")
        else:
            st.sidebar.error("Failed to configure global parameters")

    col1, col2, col3 = st.columns(3)

    with col1:
        st.subheader("Analyze by Parameter")
        group = st.selectbox("Select Group", ["Group1", "Group2", "Group3"])
        coeff_name = st.text_input("Coefficient Name")
        coeff_order = st.selectbox("Order", ["LO", "NLO", "NNLO"])
        if st.button("Retrieve Coefficient"):
            response = requests.get(
                f"{BASE_API_URL}/calculate_coefficient",
                params={"group": group, "name": coeff_name, "order": coeff_order},
            )
            if response.status_code == 200:
                data = response.json()
                st.write(f"Group: {group}")
                st.write(f"Coefficient: {coeff_name}")
                st.write(f"Value: {data['value']}")
            else:
                st.error("Failed to retrieve coefficient")

        st.subheader("Plot Variation")
        param_name = st.text_input("Parameter Name")
        param_min = st.number_input("Min Value", step=0.1)
        param_max = st.number_input("Max Value", step=0.1)
        param_steps = st.number_input("Steps", min_value=2, step=1, value=10)
        if st.button("Plot Parameter Variation"):
            response = requests.get(
                f"{BASE_API_URL}/plot_coefficients",
                params={
                    "group": group,
                    "name": coeff_name,
                    "order": coeff_order,
                    "param_name": param_name,
                    "min_value": param_min,
                    "max_value": param_max,
                    "steps": param_steps,
                },
            )
            if response.status_code == 200:
                plot_data = response.json()
                st.line_chart([val["coefficient"] for val in plot_data["values"]])
            else:
                st.error("Failed to generate plot")

    with col2:
        st.subheader("Analyze by Scale Q")
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
