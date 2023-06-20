"""
Build a Person object and have it introduce itself.

input: a string of arbitrary length, which is used to name the new person object
output: greeting printed to screen

object indroduce takes strings from the variables from class person and prints them together.
"""


class Person:
    def __init__(self, name, username, studenType, major, experience):
        self.myName = name
        self.username = username
        self.studentType = studentType
        self.major = major
        self.experience = experience

    def introduce(self):
        print('Hi there, I am {0}.'.format(self.myName))
        print('My user name is {0}.'.format(self.username))
        print('I am an {0}.'.format(self.studentType))
        print('My major is {0}.'.format(self.major))
        print('I am interested in using computers to analyse microbiology data!')


name = 'Ryan Wanless'
username = 'rwanless'
studentType = 'undergraduate'
major = 'MCD bio'
experience = 'I have zero experince with programming'
newperson = Person(name, username, studentType, major, experience)
newperson.introduce()

