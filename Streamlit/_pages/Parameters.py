import streamlit as st
import requests
import os
import matplotlib.pyplot as plt

if "wide_mode" not in st.session_state:
    st.set_page_config(layout="wide", page_title="Hyperiso", page_icon="📊")
    st.session_state["wide_mode"] = True

from Streamlit.Utils.common_elements import add_header, add_footer, apply_custom_background, apply_sidebar_style, apply_file_management_style
from Streamlit.Utils.common_elements import apply_custom_css
from Streamlit.Utils.API_utils import get_param_code_list, get_param_code_list_hist

BASE_API_URL = "http://127.0.0.1:8000"

models = {"SM", "SUSY", "THDM", "CUSTOM"}

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

if "paramtype_options" not in st.session_state:
    st.session_state.paramtype_options = None

if "hist_code_list" not in st.session_state:
    st.session_state.hist_code_list = None

if "pdg_code_list" not in st.session_state:
    st.session_state.pdg_code_list = None

def list_lha_files(directory="Assets/lha/"):
    if not os.path.exists(directory):
        return []
    return [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]

def on_new_infos():
    if st.session_state.selected_file == None:
        return
    response = requests.post(
        f"{BASE_API_URL}/parameters/set_memory_manager",
        json={"lha_file": os.path.join("Assets/lha",st.session_state.selected_file), "model": st.session_state.selected_model,
              "use_marty": st.session_state.use_marty,
        "is_spectrum": st.session_state.is_spectrum,
        "has_wilsons": st.session_state.has_wilson,
        "has_obs": st.session_state.has_obs}
    )

    if response.status_code == 200:
        st.session_state.api_status = "success"
        st.session_state.api_message = "API appelée avec succès."
    else:
        st.session_state.api_status = "error"
        st.session_state.api_message = f"Erreur lors de l'appel API : {response.status_code}"

def on_new_slha():
    selected_file = st.session_state.selected_file
    model = st.session_state.selected_model
    print(os.getcwd())
    print(os.path.join("Assets/lha", selected_file))
    if not os.path.exists("Assets/lha",selected_file):
        return
    response = requests.post(
        f"{BASE_API_URL}/parameters/set_memory_manager",
        json={"lha_file": os.path.join("Assets/lha",selected_file), "model": model, "use_marty": False,
        "is_spectrum": False,
        "has_wilsons": False,
        "has_obs": False}
    )
    print(response)
    print(response.status_code)
    if response.status_code == 200:
        st.session_state.api_status = "success"
        st.session_state.api_message = "API appelée avec succès."
    else:
        st.session_state.api_status = "error"
        st.session_state.api_message = f"Erreur lors de l'appel API : {response.status_code}"

    print(selected_file)

def on_new_model():
    print("model is: ", st.session_state.selected_model)
    selected_file = st.session_state.selected_file
    model = st.session_state.selected_model
    response = requests.post(
        f"{BASE_API_URL}/parameters/set_lha_model",
        json={"lha_file": os.path.join("Assets/lha",selected_file), "model": model, "use_marty": False,
        "is_spectrum": False,
        "has_wilsons": False,
        "has_obs": False}
    )

def app():

    apply_custom_background()
    apply_file_management_style()
    apply_custom_css()
    apply_sidebar_style()
    st.title("SLHA File Manager and Visualizer")

    st.sidebar.header("File Management")
    uploaded_file = st.sidebar.file_uploader("Upload SLHA File", type=["slha", "lha", "flha"],)
    if uploaded_file:
        response = requests.post(
            f"{BASE_API_URL}/parameters/upload",
            files={"file": (uploaded_file.name, uploaded_file.getvalue())}
        )
        if response.status_code == 200:
            st.sidebar.success("File uploaded successfully!")
        else:
            st.sidebar.error("Failed to upload file")

    if st.sidebar.button("Clean Directory"):
        lha_files = list_lha_files()
        for f in lha_files:
            os.remove(os.path.join("Assets/lha/", f))
        st.sidebar.success("Directory cleaned!")

    if not st.session_state.paramtype_options:
        response = requests.get(
                f"{BASE_API_URL}/parameters/all_blocks_list")
        data = response.json()
        print(data)
        bl = data['blocks']
        st.session_state.paramtype_options = bl

    col_left, col_center, col_right = st.columns([2.5, 2, 2])

    with col_left:
        st.subheader("SLHA File Viewer")
        lha_files = list_lha_files()
        st.session_state.api_called_for = False
        selected_file = st.selectbox("Select a SLHA File", lha_files, key= "selected_file", on_change=on_new_infos)
        selected_model = st.selectbox("Select a model", models, key = 'selected_model', on_change=on_new_infos)
        if selected_model:
            st.text(f"Current model : {selected_model}")
        if selected_file:
            file_path = os.path.join("Assets/lha/", selected_file)
            st.text(f"Contents of {selected_file}:")
            with open(file_path, "r") as f:
                st.code(f.read(), language="plaintext")

    with col_center:
        st.subheader("Single File Analysis")
        # block = st.text_input("Block Name", placeholder="Enter block name (e.g., MASS)")
        param_block = st.selectbox("Parameter Block", st.session_state.paramtype_options, on_change=get_param_code_list, key = "param_block")
        code = st.selectbox("Parameters Code", st.session_state.pdg_code_list, key = "code")
        # code = st.number_input("Parameter Code", step=1, min_value=0)
        if st.button("Get Parameter Value") and param_block and code:
            response = requests.get(f"{BASE_API_URL}/parameters/value", params={"block": param_block, "code": code})
            if response.status_code == 200:
                value = response.json().get("value")
                if value != "":
                    st.write(f"Value for block `{param_block}` and code `{code}`: {value}")
                else:
                    st.write(f"No value for block {param_block} and pdgcode {code}")
            else:
                st.error(response.json().get("detail", "Failed to retrieve parameter value"))

        st.session_state.show_pie = True
        if st.button("Show Block Distribution"):
            st.session_state.show_pie = True

        if st.session_state.show_pie:
            response_blocks = requests.get(f"{BASE_API_URL}/parameters/blocks_list", params={"param_type" : st.session_state.param_type})
            if response_blocks.status_code == 200:
                blocks = response_blocks.json().get("blocks", [])
                sizes = []
                for block in blocks:
                    response_info = requests.get(f"{BASE_API_URL}/parameters/block_info", params={"block": block, "param_type" : st.session_state.param_type})
                    if response_info.status_code == 200:
                        sizes.append(len(response_info.json().get(block, {})))
                fig, ax = plt.subplots()
                ax.pie(sizes, labels=blocks, autopct=lambda p: f'{int(p * sum(sizes) / 100)}', startangle=90)
                ax.axis("equal")
                st.pyplot(fig)
            else:
                st.error("Failed to retrieve blocks information.")
        st.write("### Deuxième graphique circulaire")
        if (not st.session_state.has_wilson) and (st.session_state.selected_model == "SM"):
            options = ["FLAVOR", "FF"]
        elif (not st.session_state.has_wilson) and (st.session_state.selected_model != "SM"):
            options = ["MODEL", "FLAVOR", "FF"]
        elif st.session_state.has_wilson and (st.session_state.selected_model == "SM"):
            options = ["FLAVOR", "FF", "WILSON"]
        else:
            options = ["MODEL", "FLAVOR", "FF", "WILSON"]

        selected_option = st.selectbox("Choisissez une option", options)

        if st.button("Générer le diagramme"):
            if selected_option:
                param_type = None
                if selected_option == "MODEL":
                    param_type = st.session_state.selected_model
                elif selected_option == "FLAVOR":
                    param_type = "FLAVOR"
                elif selected_option == "FF":
                    param_type = "FF"
                elif selected_option == "WILSON":
                    param_type = "WILSON"

                if param_type:
                    response_blocks = requests.get(f"{BASE_API_URL}/parameters/blocks_list", params={"param_type": param_type})
                    if response_blocks.status_code == 200:
                        blocks = response_blocks.json().get("blocks", [])
                        sizes = []
                        for block in blocks:
                            response_info = requests.get(f"{BASE_API_URL}/parameters/block_info", params={"block": block, "param_type": param_type})
                            if response_info.status_code == 200:
                                sizes.append(len(response_info.json().get(block, {})))

                        fig, ax = plt.subplots()
                        ax.pie(sizes, labels=blocks, autopct=lambda p: f'{int(p * sum(sizes) / 100)}', startangle=90)
                        ax.axis("equal")
                        st.pyplot(fig)
                    else:
                        st.error("Failed to retrieve blocks information for the selected option.")
                else:
                    st.warning("Veuillez sélectionner un paramètre valide.")
            else:
                st.warning("Veuillez choisir une option avant de générer le diagramme.")
    with col_right:
        st.subheader("Multi-File Analysis")

        st.session_state.show_histogram = True
        if st.button("Generate Histogram"):
            st.session_state.show_histogram = True
            
        if st.session_state.show_histogram:
            # selected_block = st.text_input("Histogram Block Name", placeholder="Enter block name (e.g., MASS)")
            # selected_code = st.number_input("Histogram Parameter Code", step=1, min_value=0)

            selected_block = st.selectbox("Parameter Block", st.session_state.paramtype_options, on_change=get_param_code_list_hist, key = "selected_block")
            selected_code = st.selectbox("Parameters Code", st.session_state.hist_code_list, key = "selected_code")
            values = []
            old_value = st.session_state.selected_file

            for lha_file in list_lha_files():
                model = "SM"
                if not selected_block or not selected_code:
                    continue
                response = requests.post(
                    f"{BASE_API_URL}/parameters/set_memory_manager",
                    json={"lha_file": os.path.join("Assets/lha",lha_file), "model": model, "use_marty": False,
                    "is_spectrum": False,
                    "has_wilsons": False,
                    "has_obs": False}
                )
                if response.status_code != 200:
                    print("Error.")
                response = requests.get(
                    f"{BASE_API_URL}/parameters/value",
                    params={"block": selected_block, "code": selected_code}
                )
                if response.status_code == 200:
                    if response.json().get("value") == "":
                        st.write(f"No value for block {selected_block} and pdgcode {selected_code}")
                        break
                    values.append(response.json().get("value"))
            if values:
                fig, ax = plt.subplots()
                ax.hist(values, bins=10, edgecolor="black")
                ax.set_title(f"Distribution of {selected_block}/{selected_code}")
                ax.set_xlabel("Value")
                ax.set_ylabel("Frequency")
                st.pyplot(fig)
            else:
                st.warning("No values available for the selected block and code.")

            on_new_infos()
    add_footer()  # Add the footer at the bottom

app()


