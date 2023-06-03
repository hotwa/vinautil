#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file               :covalent_coord
@Description:       :
@Date               :2023/4/9 17:35:40
@Author             :lyzeng
@mail               :pylyzeng@gmail.com       
@version            :1.0
'''

class CustomError(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return f"CustomError: {self.message}"
