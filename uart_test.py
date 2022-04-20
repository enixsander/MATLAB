# Terminal
# .\uart.py -s COM4 -d .\test.txt -n 512 -b 115200 -i 100
# -n - capture_num, -i - interval

# Visual Studio Code
# в файл .vscode\launch.json добавить параметры (для отладки в VS Code)
# "args": ["-s=COM4", "-d=./adc.txt", "-n=1","-f","-i=1"]

# Сначала передаём по UART 2 байта 0x27 и 0xFF, далее 512 значений (int16_t) с каждого канала АЦП (512*size(int16_t)*CHANNELS байт)
# Парсит в файл adc.txt в формате:
# -30 -33 -36 -37 -38 -40 -41 ... -36\n
# 18 21 22 23 23  ... -5\n

# Сухая проливка
# ssh test@192.168.10.168
# qwe
# cd workspace/servopack-driver
# python demo.py -t 600 -r 600

import time
import struct
import math
import os

import numpy as np
import serial
import argparse


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


class UARTPacket:

    def __init__(self, adc_buff_size):
        self.adc_buff_size = adc_buff_size
        self.adc1_buff = []
        self.adc2_buff = []
        self.adc3_buff = []
        self.adc4_buff = []
        self.save_packet_format = f"<{adc_buff_size * 2}H"

    @staticmethod
    def from_serial(serial_conn):
        p = UARTPacket(1536)
        p.update(serial_conn)
        return p

    def update(self, serial_conn):
        magic_state = 0
        magic_error_count = 0
        while magic_state != 39:
            magic_byte = serial_conn.read()
            if magic_state == 0:
                if magic_byte == b'\x27':
                    magic_state = 1
            elif magic_state == 1:
                if magic_byte == b'\xff':
                    magic_state = 39
                else:
                    magic_error_count += 1
                    magic_state = 0

        packet_opt_format = f"<{self.adc_buff_size * 4}h"

        raw_bytes = serial_conn.read(struct.calcsize(packet_opt_format))

        adc_buffs = tuple(chunks(raw_bytes, self.adc_buff_size*struct.calcsize('h')))
        adc1_buff = struct.unpack(f"{self.adc_buff_size}h", adc_buffs[0])
        adc2_buff = struct.unpack(f"{self.adc_buff_size}h", adc_buffs[1])
        adc3_buff = struct.unpack(f"{self.adc_buff_size}h", adc_buffs[2])
        adc4_buff = struct.unpack(f"{self.adc_buff_size}h", adc_buffs[3])

        for i in range(0, len(adc1_buff)):
            self.adc1_buff.append(adc1_buff[i])
            self.adc2_buff.append(adc2_buff[i])
            self.adc3_buff.append(adc3_buff[i])
            self.adc4_buff.append(adc4_buff[i])


    def to_bytes(self):
        return struct.pack(self.save_packet_format,
                           *(self.adc1_buff + self.adc2_buff))

    def to_string(self):
        adc1_str = ' '.join(str(value) for value in self.adc1_buff if value is not None)
        adc2_str = ' '.join(str(value) for value in self.adc2_buff if value is not None)
        adc3_str = ' '.join(str(value) for value in self.adc3_buff if value is not None)
        adc4_str = ' '.join(str(value) for value in self.adc4_buff if value is not None)
        data = adc1_str + '\n' + adc2_str + '\n' + adc3_str + '\n' + adc4_str + '\n' #'adc1 = [' + 
        return data

    def is_adc_12_ok(self):
        for i in self.adc1_buff:
            if i > 4096:
                return False
        return True

    def is_adc_21_ok(self):
        for i in self.adc2_buff:
            if i > 4096:
                return False
        return True

def write_to_file(data_arr, file):
    for raw_data in data_arr:
        file.write(raw_data.to_string())
    file.flush()


def write_file_header(current_time, file):
    magic_value = b"\x4d\x49\x4b\x55"
    version = 2
    current_time_ms = math.floor(current_time * 1000)
    file.write(struct.pack("<4sIQ", magic_value, version, current_time_ms))


def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument('--serial', '-s', dest='serial_path', required=True)
    argparser.add_argument('--dest', '-d', dest='destination_path', required=True)
    argparser.add_argument('--force', '-f', dest='force', action='store_true')
    argparser.add_argument('-n', dest='capture_num', required=True, type=int)
    argparser.add_argument('--baud-rate', '-b', dest='baudrate', required=False, type=int, default=115200)
    argparser.add_argument('--save-interval', '-i', dest='interval', required=False, type=int, default=100)
    args = argparser.parse_args()

    tty_path = args.serial_path
    dest_path = args.destination_path
    baudrate = args.baudrate
    capture_amount = int(args.capture_num)
    save_interval = int(args.interval)

    if os.path.exists(dest_path) and not args.force:
        raise argparse.ArgumentError(None, f'file "{dest_path}" already exists')

    with serial.Serial(tty_path, baudrate) as s, open(dest_path, 'w') as f:
        script_start_time = time.time()
        captured_packets = []

        captured_total = 0

        i = 0
        #write_file_header(script_start_time, f)
        current_packet_time = time.time()
        last_saved_time_before = current_packet_time
        while i < capture_amount:
            packet = UARTPacket.from_serial(s)
            prev_packet_time = current_packet_time
            current_packet_time = time.time()

            if not packet.is_adc_12_ok():
                print(f'cnt={i} ADC12=error')
            if not packet.is_adc_21_ok():
                print(f'cnt={i} ADC21=error')

            if 1: #packet.crc_valid():
                captured_packets.append(packet)
                i += 1
            else:
                print(f"cnt={i} CRC_ERROR")
            if i % save_interval == 0 and i > 0:
                last_saved_time = last_saved_time_before
                last_saved_time_before = time.time()
                write_to_file(captured_packets, f)
                captured_packets.clear()
                last_saved_time_after = time.time()
                capture_time = last_saved_time_before - last_saved_time
                save_time = last_saved_time_after - last_saved_time_before

                print(f'cnt={i}',
                      'avg_dt={:.3f}ms'.format(capture_time * 1000 / save_interval),
                      'save_time={:.6f}ms'.format(save_time * 1e3))

            captured_total += 1

        if len(captured_packets) > 0:
            write_to_file(captured_packets, f)

            captured_packets.clear()
        print('ok')


if __name__ == '__main__':
    main()
