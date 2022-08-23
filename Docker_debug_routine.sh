sudo docker build -f Dockerfile -t aidamri:dev .
sudo docker run --name aidamri_debug -dit --rm aidamri:dev
sudo docker cp TestData/ aidamri_debug:/aida/TestData/
sudo docker exec -w /aida/TestData/testData/testData/DTI -it aidamri_debug ls
sudo docker exec -w /aida/bin/1_PV2NIfTiConverter -it aidamri_debug python3 pv_conv2Nifti.py -i /aida/TestData/testData/testData
sudo docker stop aidamri_debug