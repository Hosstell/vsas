cd ./ncbi.d

VERSION=$(cat version.txt)
docker manifest inspect vsas/ncbi.d:$VERSION > /dev/null
if [ $? -eq 0 ]; then
  exit 0
fi

set -e

docker build -t ncbi.d .

echo "Deploying latest version of ncbi.d"
docker image tag ncbi.d vsas/ncbi.d:latest
docker image push vsas/ncbi.d:latest

VERSION=$(cat version.txt)
echo "Deploying ${VERSION} version of ncbi.d"
docker image tag ncbi.d vsas/ncbi.d:${VERSION}
docker image push vsas/ncbi.d:${VERSION}


echo "Sending notification"
cd ..
pip3 install -r ./utils/requirements.txt
python3 ./utils/telegram_notification.py -t $TELEGRAM_BOT_TOKEN -c $TELEGRAM_CHAT_ID_FOR_NOTIFICATIONS -m "$(cat ./ncbi.d/README.md)"