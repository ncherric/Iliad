#!/usr/bin/python

import yaml

def main():

    with open("config.yaml") as istream:
        ymldoc = yaml.safe_load(istream)

        # Prompt the user for a new value
        WD = input("Enter your working directory. For example, something like this -> /Path/To/Directory/Iliad/: ")
        ymldoc['workdirPath'] = WD

        # Prompt the user for a new value
        sampleFile = input("Enter your sample file location containing sample names and URLs (if you are downloading). copy this if unchanged -> config/UserSampleTable.csv: ")
        ymldoc['samplesDict'] = sampleFile

        # Prompt the user for a new value
        refSpecies = input("Enter your reference assembly species. copy this if unchanged -> homo_sapiens: ")
        ymldoc['ref']['species'] = refSpecies

        try:
            refRelease = int(input("Enter your reference assembly release. copy this if unchanged -> 104: "))
            ymldoc['ref']['release'] = refRelease
        except ValueError:
            print("Invalid input. Please enter an integer.")

        refBuild = input("Enter your reference assembly build. copy this if unchanged -> GRCh38: ")
        ymldoc['ref']['build'] = refBuild

    with open("modified.yaml", "w") as ostream:
        yaml.dump(ymldoc, ostream, default_flow_style=False, sort_keys=False)

if __name__ == "__main__":
    main()