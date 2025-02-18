import requests
import streamlit as st


BASE_API_URL = "http://127.0.0.1:8000"

def get_param_code_list():
    print(st.session_state.param_block)
    response = requests.get(
            f"{BASE_API_URL}/parameters/block_info",
            {"block" : st.session_state.param_block, "param_type" : ""})
    data = response.json()
    print(data)
    st.session_state.pdg_code_list = data[st.session_state.param_block]

def get_param_run_code_list():
    print(st.session_state.param_block_run)
    response = requests.get(
            f"{BASE_API_URL}/parameters/block_info",
            {"block" : st.session_state.param_block_run, "param_type" : ""})
    data = response.json()
    print(data)
    st.session_state.pdg_code_list = data[st.session_state.param_block_run]

def get_param_code_list_hist():
    print(st.session_state.selected_block)
    response = requests.get(
            f"{BASE_API_URL}/parameters/block_info",
            {"block" : st.session_state.selected_block, "param_type" : ""})
    data = response.json()
    print(data)
    st.session_state.hist_code_list = data[st.session_state.selected_block]
