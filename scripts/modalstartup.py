import modal

app = modal.App(name="test") 

@app.function()
def test():
	print("test")

def main():
	with app.run():
		test.remote()

