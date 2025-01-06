import streamlit as st
import requests
import os
import pandas as pd
import matplotlib.pyplot as plt

BASE_API_URL = "http://127.0.0.1:8000/parameters"

def list_lha_files(directory="DataBase/lha/"):
    """List all files in the LHA directory."""
    if not os.path.exists(directory):
        return []
    return [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]

def app():
    st.title("SLHA File Manager and Visualizer")

    # Sidebar for file management
    st.sidebar.header("File Management")
    uploaded_file = st.sidebar.file_uploader("Upload SLHA File", type=["slha", "lha"])
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

    # Main layout
    st.subheader("Choose and View SLHA Files")
    lha_files = list_lha_files()
    selected_file = st.selectbox("Select a SLHA File", lha_files)
    if selected_file:
        file_path = os.path.join("DataBase/lha/", selected_file)
        st.text(f"Contents of {selected_file}:")
        with open(file_path, "r") as f:
            st.code(f.read(), language="plaintext")

    # Blocks and parameter information
    st.subheader("Retrieve Parameter Information")
    block = st.text_input("Block Name", placeholder="Enter block name (e.g., MASS)")
    code = st.number_input("Parameter Code", step=1, min_value=0)
    if st.button("Get Parameter Value"):
        response = requests.get(f"{BASE_API_URL}/value", params={"block": block, "code": code})
        if response.status_code == 200:
            value = response.json().get("value")
            st.write(f"Value for block `{block}` and code `{code}`: {value}")
        else:
            st.error(response.json().get("detail", "Failed to retrieve parameter value"))

    # Pie chart of block distribution
    st.subheader("Block Distribution Visualization")
    if st.button("Show Block Distribution"):
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

    # Histogram of parameter values across files
    st.subheader("Parameter Distribution Across Files")
    selected_block = st.text_input("Histogram Block Name", placeholder="Enter block name (e.g., MASS)")
    selected_code = st.number_input("Histogram Parameter Code", step=1, min_value=0, key="hist_param_code")
    if st.button("Generate Histogram"):
        values = []
        for lha_file in lha_files:
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

    # Additional option: Show all blocks in a selected file
    st.subheader("Explore Blocks in a File")
    if st.button("Show Blocks in Selected File"):
        response = requests.get(f"{BASE_API_URL}/blocks_list")
        if response.status_code == 200:
            st.write("Blocks available in the file:")
            st.write(response.json().get("blocks", []))
        else:
            st.error("Failed to retrieve blocks list.")

    # Footer
    st.markdown("---")
    st.write("Developed by [Your Name] - SLHA File Manager & Visualizer")

if __name__ == "__main__":
    app()
