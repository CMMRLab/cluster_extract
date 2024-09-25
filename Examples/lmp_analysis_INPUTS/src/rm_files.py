# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
August 30th, 2022
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""


##############################
# Import Necessary Libraries #
##############################
import os


def delete(list_of_files_to_rm):
    
    
    # Find present working directory
    pwd = os.getcwd()
    
    # Remove all files in list if any file exists in list
    if list_of_files_to_rm:
        for file in list_of_files_to_rm:
            try:
                os.remove(os.path.join(pwd, file))
            except:
                pass
        
    return
