import streamlit as st
import os
def add_header():
    """Add a more polished header using Streamlit native components."""
    st.markdown(
        """
        <style>
        /* Déplacer le header en haut */
        .custom-header {
            position: absolute;
            top: 0;
            left: 0;
            width: 100%;
            background-color: #ffffff;  /* Fond blanc */
            padding: 10px 20px;
            z-index: 999; /* Toujours visible au-dessus */
            border-bottom: 2px solid #4CAF50; /* Ligne sous le header */
        }

        /* Décaler le contenu de Streamlit vers le bas pour éviter le chevauchement */
        header[data-testid="stHeader"] {
            padding-top: 80px; /* Ajuster si nécessaire */
        }
        </style>
        """,
        unsafe_allow_html=True,
    )

    # Organisation du header en colonnes pour aligner les logos et le titre
    with st.container():
        col1, col2, col3, col4 = st.columns([6, 1, 1, 1])  # Ajustement des proportions

        with col1:
            st.markdown(
                """
                <style>
                .title {
                    margin-top: -20px !important;
                    margin-bottom: 0px !important; /* Supprime l'espace sous le titre */
                }
                .subtitle {
                    font-style: italic;  /* Met en italique */
                    margin-top: -5px !important;  /* Rapproche du titre */
                    color: gray;  /* Optionnel : Mettre en gris */
                }
                
                </style>
                """,
                unsafe_allow_html=True,
            )

            st.markdown('<h3 class="title">🌈 Hyperiso</h3>', unsafe_allow_html=True)
            st.markdown('<p class="subtitle">A modern BSM flavor calculator.</p>', unsafe_allow_html=True)

            # st.write("\x1B[3mA modern BSM flavor calculator.\x1B[0m")

        with col2:
            st.image(os.path.join("Streamlit/img/cnrs.png"), width=80)

        with col3:
            st.image(os.path.join("Streamlit/img/ip2i.png"), width=80)

        with col4:
            st.image(os.path.join("Streamlit/img/mlg.jpg"), width=80)
            
    st.markdown(
        """
        <style>
        .title-line {
            border-top: 3px solid #e67f19; /* Ligne épaisse et verte */
            margin-top: 5px;
            border-radius: 5px;
            width: 98%;
        }
        </style>
        <div class="title-line"></div>
        """,
        unsafe_allow_html=True,
    )

def apply_sidebar_style(with_span = False):
    """Adjust the width of the sidebar and style its components."""
    
    # if with_span:
    #     st.markdown("""<style>
    #                 [data-testid="stSidebar"] p, 
    #     [data-testid="stSidebar"] span, 
    #     [data-testid="stSidebar"] label {
    #         color: white !important;
    #     }</style>
    #                 """, unsafe_allow_html=True,)
    # else:
    #     st.markdown("""<style>
                   
    #     }</style>
    #                 """, unsafe_allow_html=True,)
    st.markdown(
        """
        <style>
        /* Adjust the sidebar width */
        [data-testid="stSidebar"] {
            width: 350px;  /* Set the desired width */
            background-color: #213f77;
            padding: 5px; /* Add padding for better spacing */
        }
        
         [data-testid="stSidebarNav"] a {
            color: white !important;  /* Texte en blanc */
        }

        /* Modifier la couleur des liens lorsqu'on passe la souris */
        [data-testid="stSidebarNav"] a:hover {
            color: #FFFFFF !important;  /* Jaune au survol */
            border-radius: 5px;
            padding: 5px;
            color: white !important;
        }
        
        [data-testid="stSidebarNav"] span {
            color: white;
        }
        [data-testid="stSidebarNav"] {
            border-radius: 8px;
            padding: 10px;
        }

        /* Modifier la couleur du texte des options du radio */
        [data-testid="stSidebarNav"] label {
            color: white !important;  /* Texte en blanc */
            font-weight: bold !important;
        }

        /* Modifier la couleur de fond du radio sélectionné */
        [data-testid="stSidebarNav"] label[data-selected="true"] {
            color: #ffcc00 !important;  /* Jaune si sélectionné */
        }
        
        /* Style messages inside the sidebar */
        [data-testid="stSidebar"] .stAlert {
            padding-right: 30px;      /* Add padding for alerts */
            border-radius: 5px; /* Add rounded corners */
        }

        /* Optional: Customize specific alert types if desired */
        .stSuccess {
            padding-right: 400px;      /* Add padding for alerts */
            float: left;
            background-color: #d4edda; /* Light green for success */
            border-color: #c3e6cb;     /* Border color for success */
        }

        [data-testid="stSidebar"] .stError {
            background-color: #f8d7da; /* Light red for error */
            border-color: #f5c6cb;     /* Border color for error */
        }
        </style>
        """,
        unsafe_allow_html=True,
    )

def apply_custom_css():
    st.markdown(
        """
        <style>
        .stSelectbox, .stTextInput, .stNumberInput, .stButton, .success, .stAlert {
            width: 100%;
            padding: 30px;
            float: left;
        }
        </style>
        """,
        unsafe_allow_html=True
    )

def apply_custom_css_normal():
    st.markdown(
        """
        <style>
        .Selectbox, .TextInput, .NumberInput, .Button, .success, .Alert {
            width: 100%;
            padding: 30px;
            float: left;
            color: #FDF7EF;
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
        .sidebar-section {
            width: 100%;
            padding: 1px;
        }

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
            background-color: #e4e8ff !important;  /* Light gray background */
            padding: 20px;
        }

        /* Specific style for columns to keep them white */
        .stVerticalBlock {
            background-color: white;  /* White background for columns */
            padding: 10px;
            border-radius: 8px;
            # box-shadow: 0 4px 10px rgba(0, 0, 0, 0.1);  /* Add shadow for columns */
            # margin-bottom: 20px;  /* Add space between columns */
        }
        </style>
        """,
        unsafe_allow_html=True,
    )

def upper_option_streamlit_custom():
    st.markdown(
        """
        <style>
        /* Modifier la toolbar de Streamlit */
        header[data-testid="stHeader"] {
            background: rgba(0, 0, 0, 0) !important; /* Fond semi-transparent noir */
            # backdrop-filter: blur(5px); /* Flou pour un effet plus moderne */
            # height: 40px !important; /* Réduire la hauteur */
            transition: all 0.3s ease-in-out;
        }

        /* Réduire la taille du texte et icônes */
        header[data-testid="stHeader"] * {
            font-size: 12px !important;
            color: white !important; /* Texte blanc */
        }

        /* Effet au survol pour la rendre plus visible */
        header[data-testid="stHeader"]:hover {
            background: rgba(0, 0, 0, 0.7) !important;
            backdrop-filter: blur(8px);
        }
        </style>
        """,
        unsafe_allow_html=True,
    )
    
    
def vertical_space_reduction():
    st.markdown(
        """
        <style>
        /* Réduire l'espacement entre les inputs */
        div[data-testid="stVerticalBlock"] {
            margin-bottom: 0px !important;  /* Enlève l'espace en bas */
            padding-bottom: 0px !important;
        }

        /* Réduire les marges entre les widgets */
        section[data-testid="stSidebar"] div[data-testid="stVerticalBlock"] {
            margin-bottom: 5px !important;  /* Met seulement un petit espace */
        }

        /* Réduire la marge des sliders */
        div[data-baseweb="slider"] {
            margin-top: -10px !important;
            margin-bottom: -10px !important;
        }

        /* Ajuster l'espacement des labels */
        label {
            margin-bottom: -5px !important;
            margin-top: -60px !important;
            font-size: 14px !important;  /* Réduire la taille si nécessaire */
        }

        h1 {
            margin-bottom: 40px !important; /* Augmente l'espace sous les headers */
        }

        h2 {
            margin-bottom: 30px !important; /* Augmente l'espace sous les subheaders */
        }

        h3 {
            margin-bottom: 20px !important; /* Augmente l'espace sous les petits titres */
        }
    
        </style>
        """,
        unsafe_allow_html=True,
    )
    
    
def reduce_margin():
    st.markdown(
        """
        <style>
        /* Étendre la largeur de la page encore plus */
        .appview-container {
            max-width: 100% !important;  /* Étend le contenu à toute la largeur de l'écran */
            # padding-left: 10px !important;
            # padding-right: 10px !important;
        }


        /* Supprimer la marge du conteneur principal */
        .block-container {
            padding-left: 20px !important;
            padding-right: 20px !important;
            max-width: 100% !important;
        }

        .main {
            padding-left: 100px !important;
            padding-right: 100px !important;
        }
    
        /* Étendre les colonnes au maximum */
        div[data-testid="column"] {
            width: 90% !important;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )
    
def reduce_title_space():
    st.markdown(
        """
        <style>
        /* Réduire la marge en haut des titres */
        h1 {
            margin-top: -150px !important; /* Ajuste l'espace en haut */
        }
        </style>
        """,
        unsafe_allow_html=True,
    )