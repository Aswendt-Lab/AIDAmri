sudo docker build -f Dockerfile -t aidamri:dev .
sudo docker run \
-dit --rm \
--name aidamri \
--mount type=bind,source=$1,target=$2 \
aidamri:dev

sudo docker attach aidamri
