import streamlit as st
from st_pages import get_nav_from_toml
if "wide_mode" not in st.session_state:
    st.set_page_config(layout="wide", page_title="Hyperiso", page_icon="📊")
    st.session_state["wide_mode"] = True
    
from Streamlit.Utils.common_elements import apply_sidebar_style, apply_custom_background
apply_sidebar_style(True)
apply_custom_background()
# from Streamlit._pages import login, Generation, Parameters, Wilson, Observables, References, Last_update

nav = get_nav_from_toml(".streamlit/pages_sections.toml")

pg = st.navigation(nav)
pg.run()
# pg.run()


# PAGES = {
#     ":rainbow[Login]": login,
#     "LHA Generation" : Generation,
#     "Parameters": Parameters,
#     "Wilson Coefficients": Wilson,
#     "Observables": Observables,
#     "Last update": Last_update,
#     "References": References
# }

# if "authenticated" not in st.session_state:
#     st.session_state["authenticated"] = False

# apply_sidebar_style(True)
# st.sidebar.title("Navigation")
# if not st.session_state["authenticated"]:
#     st.warning("Veuillez vous connecter pour accéder à l'application.")
#     login.app()
# else:
#     apply_sidebar_style(True)
#     selection = st.sidebar.radio("Go to", list(PAGES.keys()))
#     page = PAGES[selection]
#     if hasattr(page, "app"):
#         page.app()
#     else:
#         st.error(f"La page sélectionnée n'a pas de fonction `app` : {selection}")
