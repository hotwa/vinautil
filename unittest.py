# 导入unittest模块
import unittest

# 定义一个需要测试的函数
def add(x, y):
    return x + y

# 创建一个测试类，继承自unittest.TestCase
class TestAdd(unittest.TestCase):

    # 所有的测试方法都必须以"test_"开头
    def test_add(self):
        # 使用assertEqual方法，验证add函数的结果是否正确
        self.assertEqual(add(1, 2), 3)

    def test_add_negative(self):
        # 再添加一个测试，验证负数的情况
        self.assertEqual(add(-1, -2), -3)

# 运行测试
if __name__ == '__main__':
    unittest.main()
