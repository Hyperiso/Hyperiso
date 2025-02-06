import streamlit as st
import requests
import os
import matplotlib.pyplot as plt

if "wide_mode" not in st.session_state:
    st.set_page_config(layout="wide", page_title="Hyperiso", page_icon="📊")
    st.session_state["wide_mode"] = True
    
from Streamlit.Utils.common_elements import add_header, add_footer, apply_custom_background, apply_sidebar_style, apply_file_management_style
from Streamlit.Utils.common_elements import apply_custom_css

BASE_API_URL = "http://127.0.0.1:8000/generation"

if "reference_file" not in st.session_state:
    st.session_state.reference_file = None

if "parameter_settings" not in st.session_state:
    st.session_state.parameter_settings = []

def list_lha_files(directory="DataBase/lha/"):
    if not os.path.exists(directory):
        return []
    return [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]

def check_parameter(block, number, start, stop, step):
    """Custom parameter validation."""
    if start >= stop:
        return False, "Start value must be less than stop value."
    if step <= 0:
        return False, "Step value must be greater than 0."
    return True, ""

def app():
    apply_custom_background()
    apply_file_management_style()
    apply_custom_css()
    apply_sidebar_style()
    
    st.title("SLHA File Manager and Visualizer")

    st.sidebar.header("Manage SLHA Files")
    slha_files = list_lha_files()
    selected_file = st.sidebar.selectbox("Available SLHA Files", options=slha_files)
    st.sidebar.write(f"Total files: {len(slha_files)}")
    if st.sidebar.button("Clean Directory"):
        for file in slha_files:
            os.remove(os.path.join("DataBase/lha/", file))
        st.sidebar.success("Directory cleaned successfully!")

    st.write("### Upload SLHA Template File")
    uploaded_file = st.file_uploader("Upload a SLHA template file", type=["txt", "lha", "slha", "flha"])
    if uploaded_file:
        st.session_state.reference_file = uploaded_file.read().decode("utf-8")
        st.success("Template file uploaded successfully!")

    if not st.session_state.reference_file:
        st.warning("Please upload a SLHA template file to proceed.")
        add_footer()
        return

    st.write("### Parameter Configuration")

    with st.form("add_parameter_form"):
        st.write("Add a new parameter to vary:")
        col1, col2 = st.columns(2)
        with col1:
            block = st.text_input("Block Name (e.g., MASS)", key="block_name")
            number = st.number_input("Parameter Number", min_value=1, step=1, key="parameter_number")
        with col2:
            start = st.number_input("Start Value", key="start_value")
            stop = st.number_input("Stop Value", key="stop_value")
            step = st.number_input("Step Value", min_value=0.01, step=0.01, key="step_value")
        add_button = st.form_submit_button("Add Parameter")

        if add_button:
            valid, message = check_parameter(block, number, start, stop, step)
            if not valid:
                st.error(f"Invalid parameter: {message}")
            else:
                st.session_state.parameter_settings.append({
                    "block": block,
                    "number": number,
                    "start": start,
                    "stop": stop,
                    "step": step
                })
                st.success(f"Parameter added: {block}({number}) from {start} to {stop} by {step}")

    if st.session_state.parameter_settings:
        st.write("### Current Parameter Settings")
        for idx, param in enumerate(st.session_state.parameter_settings):
            st.write(f"{idx + 1}. Block: {param['block']}, Parameter: {param['number']}, Range: {param['start']} to {param['stop']} by {param['step']}")

    if st.button("Generate SLHA Files"):
        if not st.session_state.reference_file:
            st.error("Please upload a template file before generating SLHA files.")
        elif not st.session_state.parameter_settings:
            st.error("Please add at least one parameter before generating SLHA files.")
        else:
            try:
                payload = {
                    "template": st.session_state.reference_file,
                    "parameters": st.session_state.parameter_settings
                }
                response = requests.post(f"{BASE_API_URL}/generate_slha", json=payload)
                if response.status_code == 200:
                    st.success("SLHA files generated successfully!")
                else:
                    st.error(f"Failed to generate files: {response.text}")
            except Exception as e:
                st.error(f"Error connecting to the API: {e}")

    add_footer()

app()