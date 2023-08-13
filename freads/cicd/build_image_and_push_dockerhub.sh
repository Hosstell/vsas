set -e

cd ./freads
docker build -t freads .

echo "Deploying latest version"
docker image tag freads vsas/freads:latest
docker image push vsas/freads:latest

VERSION=$(cat version.txt)
echo "Deploying ${VERSION} version"
docker image tag freads vsas/freads:${VERSION}
docker image push vsas/freads:${VERSION}


echo "Sending notification"
cd ..
pip3 install -r ./utils/requirements.txt
python3 ./utils/telegram_notification.py -t ${{ secrets.TELEGRAM_BOT_TOKEN }} -c ${{ secrets.TELEGRAM_CHAT_ID_FOR_NOTIFICATIONS }} -m "$(cat ./freads/README.md)"