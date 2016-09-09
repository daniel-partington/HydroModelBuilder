"""
Configuration Loader.

Adapted from Integrated Framework - Core Package
"""

import os
import sys
import json


class ConfigLoader():
    """
    Configuration Loader.

    Loads model configuration setup from json files
    """

    def __init__(self):
        """Configuration Loader Constructor."""
        self.model_config = {}
        sep = os.path.sep

        print Config.__file__

        try:
            self.model_config = json.load(
                open(os.path.join(os.path.dirname('__file__'), 'model_config.json')))

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
                self.model_config["working_dir"] = os.path.dirname('__file__')
            # End try

        except Exception as e:
            # If any other error occurs, default to old behaviour
            self.model_config["working_dir"] = os.path.dirname('__file__')

            print e

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

# End ConfigLoader()

# Create global variable for use
CONFIG = ConfigLoader()
