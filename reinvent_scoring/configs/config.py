import argparse
import json
import os
from pathlib import Path


DEFAULT_BASE_CONFIG_PATH = (Path(__file__).parent / 'test_config.json').resolve()

parser = argparse.ArgumentParser(description='Reinvent Scoring configuration parser')
parser.add_argument(
    '--base_config', type=str, default=DEFAULT_BASE_CONFIG_PATH,
    help='Path to basic configuration for Reinvent Scoring environment.'
)


def read_json_file(path):
    try:
        with open(path) as f:
            json_input = f.read().replace('\r', '').replace('\n', '')
        try:
            return json.loads(json_input)
        except (ValueError, KeyError, TypeError) as e:
            print(f"JSON format error in file {path}: \n {e}")
            return {}
    except FileNotFoundError:
        print(f"Configuration file not found: {path}")
        print("Using default empty configuration.")
        return {
            "DEVELOPMENT_ENVIRONMENT": False,
            "MAIN_TEST_PATH": "tmp_test_folder",
            "COMPONENT_SPECIFIC": {},
            "ENVIRONMENTAL_VARIABLES": {}
        }


args, _ = parser.parse_known_args()

reinvent_scoring_config = read_json_file(args.base_config)

# Ensure required keys exist
if "ENVIRONMENTAL_VARIABLES" not in reinvent_scoring_config:
    reinvent_scoring_config["ENVIRONMENTAL_VARIABLES"] = {}

# Set environment variables if they exist
for key, value in reinvent_scoring_config['ENVIRONMENTAL_VARIABLES'].items():
    if value:  # Only set if value is not empty
        os.environ[key] = value
