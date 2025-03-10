image="ghcr.io/varfish-org/annonars:0.41.3"
source_path="/usr/local/bin/annonars"
destination_path="/usr/local/bin/annonars"

container_id=$(docker create "$image")
docker cp "$container_id:$source_path" "$destination_path"
docker rm "$container_id"
