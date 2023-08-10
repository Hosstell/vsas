import argparse

import telebot

argParser = argparse.ArgumentParser()
argParser.add_argument("-t", "--token", type=str, required=True, help="Chat bot token")
argParser.add_argument("-c", "--chat_id", type=str, required=True, help="Chat id")
argParser.add_argument("-m", "--message", type=str, required=True, help="Message text")
args = argParser.parse_args()

TOKEN = '6572533480:AAFoasjyoQOEobAkOywngiQWVYUEhApEZkc'
CHAT_ID = '-1001753512967'


TOKEN = args.token
CHAT_ID = args.chat_id
MESSAGE_TEXT = args.message

bot = telebot.TeleBot(token=TOKEN)
bot.send_message(CHAT_ID, MESSAGE_TEXT, parse_mode='Markdown')
