from os import PathLike
import yaml


def read_config(file_path: str | PathLike) -> dict:
    """
    Reads a YAML configuration file and returns the data as a dictionary.

    Args:
    file_path (str): The path to the YAML file to be read.

    Returns:
    dict: The parsed YAML data.
    """
    try:
        with open(file_path, "r") as file:
            config = yaml.safe_load(file)
            return config
    except FileNotFoundError:
        raise FileNotFoundError(f"The file {file_path} was not found.")
    except yaml.YAMLError as exc:
        print(f"Error while parsing YAML file: {exc}")
    except Exception as e:
        print(f"An error occurred: {e}")
