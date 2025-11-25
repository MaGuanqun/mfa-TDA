set -e   # stop if any script fails

echo "Running first script..."
cd ./CoordNet/Code/
bash coordnet_server.sh

echo "Running second script..."

cd ../..
bash critical_point_tracking_INR.sh

echo "All scripts finished."