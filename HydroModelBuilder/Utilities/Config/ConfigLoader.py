"""
Configuration Loader.
"""

import os
import sys
from jsmin import jsmin
import json
import subprocess

class ConfigLoader(object):
    """
    Configuration Loader.

    Loads model configuration setup from json files
    """

    def __init__(self, config_file=None):
        """Configuration Loader Constructor."""
        self.model_config = {}
        sep = os.path.sep

        if config_file is not None:
            here = os.path.dirname(config_file)
            try:
                with open(config_file) as _invalidJSON:
                    temp = jsmin(_invalidJSON.read());
                # End with
            except Exception:
                sys.exit("Error reading config file. Could be due to incorrect JSON formatting")
            # End try

            self.model_config = json.loads(temp)
            return
        else:
            # Maintain previous behaviour
            # Open config file from wherever current directory is
            # here = os.path.dirname(__file__)
            here = "" # os.path.dirname(__file__)
            config_file = 'model_config.json'
        # End if

        try:
            with open(config_file) as _invalidJSON:
                temp = jsmin(_invalidJSON.read());
                self.model_config = json.loads(temp)
            # End with

            # Check that working directory exists and ends with slash
            try:

                if str(self.model_config["working_dir"]) != "":

                    if not self.model_config["working_dir"].endswith(sep):
                        # Add missing ending slash if not provided
                        self.model_config["working_dir"] = self.model_config["working_dir"] + sep
                    # End if

                    print "Working Directory set to: " + str(self.model_config["working_dir"])

                    # Change to defined working directory
                    os.chdir(self.model_config["working_dir"])
                # End if

            except KeyError:
                sys.exit("ERROR: Invalid working directory or working directory is not set")
            except OSError:
                print("WARNING: Invalid working directory, defaulting to current directory")

                # If any error occurs, default to old behaviour
                self.model_config["working_dir"] = here
            # End try

        except Exception as e:
            # If any other error occurs, default to old behaviour
            self.model_config["working_dir"] = here

            # Hide exception msg if this has resorted to default behaviour
            if config_file == here:
                print e
            # End if

        # End try

        if 'relative_paths' in self.model_config:
            # Build paths
            for name, path in self.model_config['relative_paths'].iteritems():

                self.model_config[name] = self.model_config["working_dir"] + path

                if not self.model_config[name].endswith(sep):
                    self.model_config[name] = self.model_config[name] + sep
                # End if
            # End for
        # End if

    # End __init__()

    def set_environment(self, project_name, user_env_name=None):
        """
        Set intended project environment. Name must match one of those defined in the config file.

        In addition to changing configurations, this also changes directories to the project folder
        as given in the config file.

        :param user_env_name: str, name of the user environment.
                              If not given (i.e. set to None), uses the logged in username.
                              Defaults to None.
        """
        if user_env_name is None:
            user_env_name = subprocess.check_output("whoami").strip()

        self.settings = self.model_config[project_name]["environment"][user_env_name]
        self.model_config = self.model_config[project_name]["environment"][user_env_name]

        for v, var in self.settings.iteritems():
            setattr(self, v, var)
        # End for

        proj_folder = self.settings.get("project_folder", "")

        if len(proj_folder) > 0:
            os.chdir(proj_folder)
        # End if

        return self
    # End set_environment()

    def get_setting(self, params):
        """
        Retrieve value from settings.

        :param params: List, list of parameters to drill down
                       e.g. ["models", "Climate", "path"] will get value(s) for
                            CONFIG.settings["models"]["Climate"]["path"]

        :returns: object, value for specified parameter
        """
        temp = self.settings
        for p in params:
            temp = temp.get(p, {})
        # End for
        return temp
    # End get_setting()

# End ConfigLoader()

# Create global variable for use
CONFIG = ConfigLoader()
