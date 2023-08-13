cd ./freads

VERSION=$(cat version.txt)
docker manifest inspect vsas/freads:$VERSION > /dev/null
if [ $? -eq 0 ]; then
  exit 0
fi

set -e

docker build -t freads .

echo "Deploying latest version"
docker image tag freads vsas/freads:latest
docker image push vsas/freads:latest

echo "Deploying ${VERSION} version"
docker image tag freads vsas/freads:${VERSION}
docker image push vsas/freads:${VERSION}

echo "Sending notification"
cd ..
pip3 install -r ./utils/requirements.txt
python3 ./utils/telegram_notification.py -t $TELEGRAM_BOT_TOKEN -c $TELEGRAM_CHAT_ID_FOR_NOTIFICATIONS -m "$(cat ./freads/README.md)"