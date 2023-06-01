# PyProtif Plugin

PyProtif is a plugin for the PyMOL molecular visualization program; it automatically retrieves protein structural and/or functional motifs from different databases (PDBSum and Genomenet respectively) and integrates them in PyMOL for visualization and analyses. Through an expendable menu and a user-friendly interface, the plugin grants the users the ability to study simultaneously multiple proteins and to select and manipulate each motif separately.

## Installation instructions:

* Download the source code by clicking on the green 'Code' button, then 'Download'.
* Open PyMOL, and in the menu bar click Plugin -> Plugin Manager
* Click the tab 'Install new plugin' and then the button 'Choose file', select the source code you just downloaded. 
* Accept the default plugin installation directory or change to your preference

You should now have the entry *`PLUGIN NAME`* in your 'Plugin' menu then 'Legacy Plugins'. 

Python routines were added at the top of the plugin source code to install the libraries that are not found in the default setup of Python to avoid errors during the initialization of the plugin.
