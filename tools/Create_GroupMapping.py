import os
import tkinter as tk
from tkinter import filedialog

# Function to filter and list the subfolders based on the criteria
def filter_and_list_subfolders(main_folder, search_string):
    subfolder_list = []

    # Recursively iterate over the subfolders in the main folder
    for root, dirs, files in os.walk(main_folder):
        for folder_name in dirs:
            folder_path = os.path.join(root, folder_name)

            # Check if the search string is present in the folder name
            if search_string in folder_name:
                subfolder_list.append(folder_name)

    return subfolder_list

# Function to handle the button click event
def generate_list():
    # Get the selected main folder and search string
    main_folder = main_folder_var.get()
    search_string = search_string_var.get()

    # Filter and list the subfolders based on the criteria
    subfolder_list = filter_and_list_subfolders(main_folder, search_string)

    # Save the list of subfolder names to a text file
    list_file_path = filedialog.asksaveasfilename(defaultextension=".txt")
    with open(list_file_path, 'w') as list_file:
        list_file.write('\n'.join(subfolder_list))

    tk.messagebox.showinfo("List Generation Complete", "Subfolder names have been listed and saved!")

# Create the main window
window = tk.Tk()
window.title("Subfolder List Generator")

# Create variables to store the selected main folder and search string
main_folder_var = tk.StringVar()
search_string_var = tk.StringVar()

# Function to handle the "Browse" button click event for selecting the main folder
def select_main_folder():
    main_folder = filedialog.askdirectory()
    main_folder_var.set(main_folder)

# Create labels and entry fields for main folder selection and search string input
tk.Label(window, text="Main Folder:").grid(row=0, column=0)
tk.Entry(window, textvariable=main_folder_var, width=40).grid(row=0, column=1)
tk.Button(window, text="Browse", command=select_main_folder).grid(row=0, column=2)

tk.Label(window, text="Search String:").grid(row=1, column=0)
tk.Entry(window, textvariable=search_string_var, width=20).grid(row=1, column=1)

# Create a button to generate the list
tk.Button(window, text="Generate List", command=generate_list).grid(row=2, column=1)

# Start the main event loop
window.mainloop()
