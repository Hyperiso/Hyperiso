import streamlit as st
import requests
import os
import matplotlib.pyplot as plt

# Import the common elements
from Streamlit.Utils.common_elements import add_header, add_footer, apply_custom_background, apply_sidebar_style, apply_file_management_style

BASE_API_URL = "http://127.0.0.1:8000/parameters"


def list_lha_files(directory="DataBase/lha/"):
    if not os.path.exists(directory):
        return []
    return [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]

def app():

    apply_custom_background()  # Apply global styling
    add_header()          # Add the header
    apply_file_management_style()
    st.title("SLHA File Manager and Visualizer")

    st.sidebar.header("File Management")
    uploaded_file = st.sidebar.file_uploader("Upload SLHA File", type=["slha", "lha"],)
    if uploaded_file:
        response = requests.post(
            f"{BASE_API_URL}/upload",
            files={"file": (uploaded_file.name, uploaded_file.getvalue())}
        )
        if response.status_code == 200:
            st.sidebar.success("File uploaded successfully!")
        else:
            st.sidebar.error("Failed to upload file")

    if st.sidebar.button("Clean Directory"):
        lha_files = list_lha_files()
        for f in lha_files:
            os.remove(os.path.join("DataBase/lha/", f))
        st.sidebar.success("Directory cleaned!")

    col_left, col_center, col_right = st.columns([2.5, 2, 2])

    with col_left:
        st.subheader("SLHA File Viewer")
        lha_files = list_lha_files()
        selected_file = st.selectbox("Select a SLHA File", lha_files)
        if selected_file:
            file_path = os.path.join("DataBase/lha/", selected_file)
            st.text(f"Contents of {selected_file}:")
            with open(file_path, "r") as f:
                st.code(f.read(), language="plaintext")

    with col_center:
        st.subheader("Single File Analysis")
        block = st.text_input("Block Name", placeholder="Enter block name (e.g., MASS)")
        code = st.number_input("Parameter Code", step=1, min_value=0)
        if st.button("Get Parameter Value"):
            response = requests.get(f"{BASE_API_URL}/value", params={"block": block, "code": code})
            if response.status_code == 200:
                value = response.json().get("value")
                st.write(f"Value for block `{block}` and code `{code}`: {value}")
            else:
                st.error(response.json().get("detail", "Failed to retrieve parameter value"))

        st.session_state.show_pie = True
        if st.button("Show Block Distribution"):
            st.session_state.show_pie = True

        if st.session_state.show_pie:
            response_blocks = requests.get(f"{BASE_API_URL}/blocks_list")
            if response_blocks.status_code == 200:
                blocks = response_blocks.json().get("blocks", [])
                sizes = []
                for block in blocks:
                    response_info = requests.get(f"{BASE_API_URL}/block_info", params={"block": block})
                    if response_info.status_code == 200:
                        sizes.append(len(response_info.json().get(block, {})))
                fig, ax = plt.subplots()
                ax.pie(sizes, labels=blocks, autopct='%1.1f%%', startangle=90)
                ax.axis("equal")
                st.pyplot(fig)
            else:
                st.error("Failed to retrieve blocks information.")

    with col_right:
        st.subheader("Multi-File Analysis")

        st.session_state.show_histogram = True
        if st.button("Generate Histogram"):
            st.session_state.show_histogram = True
            
        if st.session_state.show_histogram:
            selected_block = st.text_input("Histogram Block Name", placeholder="Enter block name (e.g., MASS)")
            selected_code = st.number_input("Histogram Parameter Code", step=1, min_value=0)
            values = []
            for lha_file in list_lha_files():
                response = requests.get(
                    f"{BASE_API_URL}/value",
                    params={"block": selected_block, "code": selected_code}
                )
                if response.status_code == 200:
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

    add_footer()  # Add the footer at the bottom

if __name__ == "__main__":
    app()


