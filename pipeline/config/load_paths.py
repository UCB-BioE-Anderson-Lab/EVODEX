import yaml

def load_paths(config_path: str) -> dict:
    with open(config_path, 'r') as file:
        paths = yaml.safe_load(file)
    return paths
