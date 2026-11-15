The `profiles/` directory can contain any number of subdirectories, each containing a `config.yaml` file with a [profile specific to a workflow](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles):

`profiles/<specific_profile_name>/config.yaml`

The profile `profiles/default/config.yaml` will automatically be used by Snakemake whenever you do not provide a profile specific to a workflow via `--workflow-profile`.
This means that any resources or other (command line) arguments specified there, will implicitly be used when running this workflow.
Thus, as a workflow author, only put configurations there that you expect to work in most environments, but which the users might want to tweak.
For resource settings specific to a rule, preferably provide generally applicable settings in the rule definition, if necessary through a [dynamic resource](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#dynamic-resources) specification. Users can always override those settings in a profile if needed.

For any more specific profiles, use separate and clearly named subdirectories.
For example, use `profiles/slurm/config.yaml` for a profile specific to Slurm, or something like `profiles/slurm_uni_xyz/config.yaml` for a particular institutional Slurm compute cluster.

It is also good practice to add clear documentation comments for each entry in a (workflow) profile.
This should explain the respective entry, indicate what kind of values can be used and why a particular value or setting were chosen.
To this end, it is often helpful to provide links to relevant documentation pages, either from snakemake, a snakemake plugin or a specific cluster environment.

In general, we welcome pull requests that add such a profile for your specific compute environment to third party workflows you use.
