docker build -f dockerfile -t eagleshot/switch_finder_image:latest  -t eagleshot/switch_finder_image:1.0 .
docker push eagleshot/switch_finder_image:latest
docker push eagleshot/switch_finder_image:1.0
