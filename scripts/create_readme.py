import argparse


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Fill in the template with provided variables."
    )
    parser.add_argument("template_file", type=str, help="Path to the template file.")
    parser.add_argument(
        "output_file",
        type=str,
        help="Path to the output file where the filled template will be saved.",
    )
    parser.add_argument(
        "--project_title", type=str, required=True, help="The title of the project."
    )
    parser.add_argument(
        "--project_organization",
        type=str,
        required=True,
        help="The GitHub organization of the project.",
    )
    parser.add_argument(
        "--project_repo",
        type=str,
        required=True,
        help="The GitHub repository of the project.",
    )
    parser.add_argument(
        "--project_version", type=str, required=True, help="The version of the project."
    )
    parser.add_argument(
        "--project_documentation_url",
        type=str,
        required=True,
        help="The URL for the project documentation.",
    )
    return parser.parse_args()


def read_template(file_path):
    """Read the template file and return its content."""
    with open(file_path, "r") as file:
        return file.read()


def write_output(file_path, content):
    """Write the filled template to an output file."""
    with open(file_path, "w") as file:
        file.write(content)


def main():
    args = parse_arguments()

    # Read the template content
    template_content = read_template(args.template_file)

    # Prepare the replacement dictionary
    variables = {
        "project_title": args.project_title,
        "project_organization": args.project_organization,
        "project_repo": args.project_repo,
        "project_version": args.project_version,
        "project_documentation_url": args.project_documentation_url,
    }

    # Replace placeholders in the template
    filled_content = template_content.format(**variables)

    # Write the filled template to the output file
    write_output(args.output_file, filled_content)


if __name__ == "__main__":
    main()
