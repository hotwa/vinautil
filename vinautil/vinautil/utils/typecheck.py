#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file        :typecheck.py
@Description:       :
@Date     :2022/11/11 15:17:33
@Author      :hotwa
@version      :1.0
'''
# Descriptor for a type-checked attribute
class Typed:
    def __init__(self, name, expected_type):
        self.name = name
        self.expected_type = expected_type

    def __get__(self, instance, cls):
        if not instance:
            return self
        else:
            return instance.__dict__[self.name]

    def __set__(self, instance, value):
        if not isinstance(value, self.expected_type):
            raise TypeError(f'Expected {str(self.expected_type)}, not {type(value)}')
        instance.__dict__[self.name] = value

    def __delete__(self, instance):
        del instance.__dict__[self.name]


# Class decorator that applies it to selected attributes
def typeassert(**kwargs):
    def decorate(cls):
        for name, expected_type in kwargs.items():
            # Attach a Typed descriptor to the class
            setattr(cls, name, Typed(name, expected_type))
        return cls
    return decorate


@typeassert(course_name=str, total_class=int, score=float)
class Stock:
    def __init__(self, course_name, total_class, score):
        self.course_name = course_name
        self.total_class = total_class
        self.score = score


# base of use
# define the "Floater"
class Floater:
    def __init__(self, name):
        self.name = name

    def __get__(self, instance, cls):
        if not instance:
            return self
        else:
            return instance.__dict__[self.name]

    def __set__(self, instance, value):
        if not isinstance(value, float):
            raise TypeError("Expected an float object")
        instance.__dict__[self.name] = value

    def __delete__(self, instance):
        del instance.__dict__[self.name]


# Descriptor attribute for an integer type-checked attribute
class Integer:
    def __init__(self, name):
        self.name = name

    def __get__(self, instance, cls):
        if not instance:
            return self
        else:
            return instance.__dict__[self.name]

    def __set__(self, instance, value):
        if not isinstance(value, int):
            raise TypeError('Expected an int object')
        instance.__dict__[self.name] = value


class Point:
    x = Floater('x')
    y = Floater('y')

    def __init__(self, x, y):
        self.x, self.y = x, y