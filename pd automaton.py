# -*- coding: utf-8 -*-
"""
Created on Wed May 15 11:59:03 2024

@author: HUFFMP
"""
# import win32com.client as wcom
# from win32com.client import combrowse


# def start_proteome_discoverer_search():

#     # Create a connection to Proteome Discoverer
#     pd_app = wcom.Dispatch("Thermo.ProteomeDiscoverer.ConsensusWorkflow")

#     # # Load an existing workflow
#     # workflow_path = r"C:\Path\To\Your\Workflow.pdWorkflow"
#     # workflow = pd_app.Workflows.LoadWorkflow(workflow_path)

#     # # Start the search
#     # workflow.Execute()

#     # print("Search started successfully!")


# if __name__ == "__main__":
#     start_proteome_discoverer_search()

from win32com.client import combrowse
combrowse.main()