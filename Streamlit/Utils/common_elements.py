import streamlit as st
import os
def add_header():
    """Add a more polished header using Streamlit native components."""
    with st.container():
        col1, col2, col3, col4 = st.columns([5, 1, 1, 1])  # Split header into two parts

        with col1:
            st.markdown(
                """
                <style>
                .header-title {
                    font-size: 28px;
                    font-weight: bold;
                    color: #4CAF50;
                    margin-bottom: 0;
                }
                </style>
                <div class="header-title">Hyperiso</div>
                """,
                unsafe_allow_html=True,
            )
            st.write("A software for modern BSM Observable calculation. Much fun.")

        with col2:
            # Display logos in a row
            st.image(os.getcwd() + "/Streamlit/img/cnrs.png",
                width=100,
                use_container_width="auto",
            )
        with col3:
            st.image(os.getcwd() +"/Streamlit/img/ip2i.png",
                width=100,
                use_container_width="auto",
            )

        with col4:
            st.image(os.getcwd() + "/Streamlit/img/mlg.jpg"
                ,
                width=100,
                use_container_width="auto",
            )

def apply_sidebar_style():
    """Adjust the width of the sidebar to fit larger widgets."""
    st.markdown(
        """
        <style>
        /* Adjust the sidebar width */
        [data-testid="stSidebar"] {
            width: 350px;  /* Set the desired width */
            padding: 10px; /* Add padding for better spacing */
        }
        </style>
        """,
        unsafe_allow_html=True,
    )

def apply_custom_css():
    st.markdown(
        """
        <style>
        /* Restreindre la largeur des selectbox dans une colonne */
        .stSelectbox, .stTextInput, .stNumberInput {
            width: 100%; /* Limite la largeur à 100% de la colonne */
            padding: 30px;
            float: left;
        }
        </style>
        """,
        unsafe_allow_html=True
    )

def apply_file_management_style():
    """Style adjustments specifically for the File Management section."""
    st.markdown(
        """
        <style>
        /* Adjust the width of the File Management section */
        .sidebar-section {
            width: 100%;  /* Ensure the section takes the full width of the sidebar */
            padding: 1px;
        }

        /* Adjust the file uploader widget to align properly */
        [data-testid="stFileUploader"] {
            float: left;
            padding: 0;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )

def add_footer():
    """Add a styled footer with authors, email, and license."""
    st.markdown(
        """
        <style>
        .footer {
            text-align: center;
            font-size: 14px;
            color: #777;
            margin-top: 20px;
            padding-top: 10px;
            border-top: 1px solid #ddd;
        }
        </style>
        <div class="footer">
            <p>Developed by Théo Reymermier, Niels Fardeau, Nazila Mahmoudi. Contact us at <a href="t.reymermier@ip2i.in2p3.fr">t.reymermier@ip2i.in2p3.fr</a></p>
            <p>&copy; 2025 Hyperiso. Licensed under the MIT License.</p>
        </div>
        """,
        unsafe_allow_html=True,
    )

def apply_custom_background():
    """Apply a custom background color to the app with white columns."""
    st.markdown(
        """
        <style>
        /* Global background for the entire app */
        .stApp {
            background-color: #f2f2f2;  /* Light gray background */
            padding: 20px;
        }

        /* Specific style for columns to keep them white */
        .stVerticalBlock {
            background-color: white;  /* White background for columns */
            padding: 15px;
            border-radius: 8px;
            box-shadow: 0 4px 10px rgba(0, 0, 0, 0.1);  /* Add shadow for columns */
            margin-bottom: 20px;  /* Add space between columns */
        }
        </style>
        """,
        unsafe_allow_html=True,
    )
