import os

def del_file(root_dir=".",file='logo.png'):
    """
    @Function: Delete all bmp image file in root_dir and its subdirectory 
    @root_dir: The target directory 
    """
    file_list = os.listdir(root_dir)
    for f in file_list:
        file_path = os.path.join(root_dir, f)
        if os.path.isfile(file_path):
            if f.endswith(file):
                os.remove(file_path)
                print (" File removed! " + file_path)
        elif os.path.isdir(file_path):
            del_file(file_path)
