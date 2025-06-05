import os
import shutil

def setup_directories(base_dir, paths):
    images_dir = os.path.join(base_dir, 'images')
    data_dir = os.path.join(base_dir, 'data')
    pages_dir = os.path.join(base_dir, 'pages')

    os.makedirs(images_dir, exist_ok=True)
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(pages_dir, exist_ok=True)

    # Copy necessary files to the website directory
    if 'template_dir' in paths:
        os.makedirs(os.path.join(base_dir, 'templates'), exist_ok=True)
        for template_file in os.listdir(paths['template_dir']):
            src = os.path.join(paths['template_dir'], template_file)
            dst = os.path.join(base_dir, 'templates', template_file)
            if os.path.isfile(src):
                shutil.copy(src, dst)
    if 'static_dir' in paths:
        os.makedirs(os.path.join(base_dir, 'static'), exist_ok=True)
        for static_file in os.listdir(paths['static_dir']):
            src = os.path.join(paths['static_dir'], static_file)
            dst = os.path.join(base_dir, 'static', static_file)
            if os.path.isfile(src):
                shutil.copy(src, dst)
