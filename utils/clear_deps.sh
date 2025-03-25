# Clears unwanted dependencies. Should be run from root directory as:
# source utils/clear_deps.sh

# Clear unwanted dependencies from docs/Project.Toml

sed -i '/QGDipoles =/d' docs/Project.toml
sed -i '/GeophysicalFlows =/d' docs/Project.toml
sed -i '/JuliaFormatter =/d' docs/Project.toml

# Clear unwanted dependencies from Project.Toml

sed -i '/QGDipoles =/d' Project.toml
sed -i '/Plots =/d' Project.toml
sed -i '/GeophysicalFlows =/d' Project.toml
sed -i '/Documenter =/d' Project.toml
sed -i '/JuliaFormatter =/d' Project.toml
