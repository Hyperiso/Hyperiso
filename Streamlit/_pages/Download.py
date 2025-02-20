import streamlit as st
import requests
from Streamlit.Utils.common_elements import add_header, add_footer, apply_custom_background, apply_sidebar_style, apply_file_management_style
from Streamlit.Utils.common_elements import apply_custom_css
import os

if "wide_mode" not in st.session_state:
    st.set_page_config(layout="wide", page_title="Hyperiso", page_icon="📊")
    st.session_state["wide_mode"] = True
    
BASE_API_URL = "http://127.0.0.1:8000"

def download_request():
    response = requests.get(f"{BASE_API_URL}/download/download", stream=True)
    
    if response.status_code == 200:
        file_path = "hyperiso_latest.zip"
        with open(file_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        
        with open(file_path, "rb") as f:
            st.success("Download successfully initiated!")
            st.download_button(label="📥 Click to Download Hyperiso latest",
                               data=f,
                               file_name="hyperiso_latest.zip",
                               mime="application/zip")
        
        os.remove(file_path)
    else:
        st.error("Error while initiating the download.")

    

def app():
    apply_custom_background()
    apply_file_management_style()
    apply_custom_css()
    
    with st.sidebar:
        apply_sidebar_style()
        st.header("Download")
        if st.button("📥 Download a version of Hyperiso ?"):
            download_request()
        st.markdown("---")
        st.subheader("Resources")
        st.markdown("[📖 Manual](https://en.wikipedia.org/wiki/Rabbit)")
        st.markdown("[📑 GitHub Documentation](https://hyperiso.github.io/Hyperiso/)")
    
    st.title("Welcome to Hyperiso")
    st.markdown(
        """
        ## Advanced Computational Framework for flavour physics calculation 🔬
        Hyperiso is a flavor observable caclulator, integrating **Marty**, **Softsusy**, and **2HDMC** 
        for an in-depth exploration of extended Standard Model scenarios.
        
        ### 🔹 Key Features:
        - 📂 **Automated Calculations**: Calculation of wilson coefficients (B coefficients) and observables.
        - 🔍 **Parameter generation**: Softsusy and 2HDMC for Spectrum generation.
        - ⚡ **Advanced Statistics**: Chi2 calculation for combined observables analysis.
        - 🚀 **Seamless Marty Integration**: Automatic calculation of Wilson coefficients in new models.
        
        Hyperiso enables API-driven request processing to streamline and automate these tasks in the background.
        The download button allows access to the latest compiled version of the software suite.
        
        **Ready to deepen your research? Download Hyperiso now!**
        """
    )

    add_footer()

app()