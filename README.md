# PyProtif Plugin

## Introduction

PyProtif is a plugin for the PyMOL molecular visualization program; it automatically retrieves protein structural and/or functional motifs from different databases (PDBsum and GenomeNet respectively) and integrates them in PyMOL for visualization and analyses. Through an expendable menu and a user-friendly interface, the plugin grants the users the ability to study simultaneously multiple proteins and to select and manipulate each motif separately.

---

## Installation instructions

PyProtif plugin is fully written in Python programming language (version 3.x). Many routines were added at the beginning of the code to manage its dependencies on the local machine (install the libraries that are not part of the default Python setup) before integrating it into PyMOL.

To install the PyProtif plugin, follow these steps.

1. Make sure you have PyMOL installed on your system. PyProtif requires PyMOL 2.X.
2. The plugin utilizes the web scraping capabilities of Selenium libraries to retrieve the various motifs. Please ensure that the most recent stable version of Google Chrome Browser is installed on your system for optimal functionality.
3. Download the source code by clicking on the green 'Code' button, then 'Download'.
4. Open PyMOL, and in the menu bar click Plugin -> Plugin Manager
5. Click the tab 'Install new plugin' and then the button 'Choose file', select the source code you just downloaded.
6. Accept the default plugin installation directory or change to your preference.
7. Restart PyMOL

You should now have the entry `'PyProtif'` in the `'Legacy Plugins'` under the `'Plugin'` menu in PyMOL. However, make sure that the option to load the plugin on startup is selected.

### For Linux users

Since the plugin is fully written in Python 3.X and relies on the `Selenium` library and `chromedriver` for web scraping, it is imperative to have Python 3.X and the latest stable version of the Google Chrome browser installed on your system for the plugin to function properly.

To proceed, open a terminal by pressing `Ctrl + Alt + T`. Update your local system's repository list and install the latest version of Python by entering the following commands:

```
sudo apt update
sudo apt install python3
sudo apt install python3-pip
```

All the above commands must be executed with administrator permissions (sudo).

Python 3.X must be set as the default python version for your Linux system. This can be done by creating aliases for Python 3.X and PIP (the package manager for Python) which is needed for the automatic installation of any missing Python libraries required for the plugin to function correctly.

To ensure that the below commands are always available, you can add them to the `.bashrc` file. After making these changes, you can either source the file (`source ~/.bashrc`) or restart the terminal session for these newly applied configurations to take effect.

```
alias python=python3
alias pip=pip3
```

If users are still struggling with the dependencies due to local configurations and system settings that are conflicting with the proper installation/functioning of the plugin, we recommend them to install these via a terminal by typing:

```
python3 -m pip install selenium
python3 -m pip install chromedriver-autoinstaller
```

All others dependencies are default Python libraries.

At this point, users can try to reinstall the plugin using the PyMOL plugin manager, and this should result now in a successful installation.

## Usage

For an illustrative guide on how to use the PyProtif plugin, please refer to the `PyProtif_tutorial.md` file.

## Contact

If you have any questions, suggestions, bug reports or feedback, please contact us through the GitHub Discussions forum of the project.
