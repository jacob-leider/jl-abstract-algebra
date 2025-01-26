import math
import secrets

class DiffieHellmanClient:
  def __init__(self, **kwargs):
    # Bytes per block.
    self.block_size_ = 128
    if "block_size" in kwargs:
      self.block_size_ = kwargs["block_size"]
    # Mod (`p`: public).
    self.p_ = None
    if "public_prime" in kwargs:
      self.p_ = kwargs["public_prime"]
    # Base (`g`: public).
    self.g_ = None
    if "public_base" in kwargs:
      self.g_ = kwargs["public_base"]
    # Exponent (`a`: private) and my shared factor (`A` = `g`^`a`: public).
    self.private_exp_ = None
    self.A_ = None
    if kwargs.get("auto_set_private_exp", False):
      self.set_private_exp()
      self.compute_g_pow_a()
    # Your shared factor (`B` = `g`^`b`: public).
    self.B_ = None
    # Shared public key ((`g`^`a`)^`b`: public).
    self.shared_public_key_ = None

  def set_public_prime(self, p: int):
    if math.ceil(math.log2(p) / 8) < self.block_size_:
      raise Exception("prime size (in bytes) must be at least as large as block size.")
    self.p_ = p

  def set_public_base(self, g: int):
    self.g_ = g

  def set_private_exp(self):
    if self.p_ == None:
      raise Exception("mod (p_) not set")
    self.private_exp_ = secrets.randbelow(self.p_ - 1) + 1

  def compute_g_pow_a(self):
    if self.p_ == None:
      raise Exception("mod (p_) not set")
    if self.g_ == None:
      raise Exception("base (g_) not set")
    if self.private_exp_ == None:
      self.set_private_exp()
    self.A_ = pow(self.g_, self.private_exp_, self.p_)

  def set_shared_public_key(self, B: int):
    if self.p_ == None:
      raise Exception("mod (p_) not set")
    if self.g_ == None:
      raise Exception("base (g_) not set")
    if self.private_exp_ == None:
      self.set_private_exp()
    if self.A_ == None:
      self.compute_shared_public_key_factor()
    self.B_ = B
    self.shared_public_key_ = pow(B, self.private_exp_, self.p_)

  def _decrypt_block(self, msg_block: bytes, c1_pow_a_inv: int) -> bytes:
    c2_int = int.from_bytes(msg_block)
    decrypted_block_int = (c1_pow_a_inv * c2_int) % self.p_
    decrypted_block_bytes = decrypted_block_int.to_bytes(self.block_size_, byteorder='big')
    return decrypted_block_bytes

  def _encrypt_block(self, encrypted_block: bytes, k: int) -> bytes:
    encrypted_block_int = int.from_bytes(encrypted_block)
    c2_int = (encrypted_block_int * pow(self.B_, k, self.p_)) % self.p_
    c2_bytes = c2_int.to_bytes(self.block_size_, byteorder='big')
    return c2_bytes

  def encrypt_msg(self, msg: str | bytes) -> bytes:
    # Error handling.
    if self.A_ == None:
      raise Exception("source's factor (A_) not set")
    if self.B_ == None:
      raise Exception("target's factor (B_) not set")

    # Preprocess into byte string of valid size.
    if msg.__class__ not in (str, bytes):
      raise Exception("msg must be str or bytes")
    if msg.__class__ == str:
      msg = bytes(msg, encoding='utf-8') # Always UTF 8
    p_size_bytes = math.ceil(math.log2(p) / 8)
    num_encrypted_blocks = math.ceil(len(msg) / self.block_size_)
    num_encrypted_bytes = num_encrypted_blocks * self.block_size_
    encrypted_msg = bytearray(num_encrypted_bytes)

    # Generate temporary `k`.
    k = secrets.randbelow(self.p_)

    # Encrypt.
    c1_int = pow(self.g_, k, self.p_)
    for i in range(0, num_encrypted_blocks):
      encrypted_slice = slice(self.block_size_ * i, self.block_size_ * (i + 1))
      decrypted_block_bytes = msg[encrypted_slice]
      encrypted_msg[encrypted_slice] = self._encrypt_block(decrypted_block_bytes, k)

    c1_bytes = c1_int.to_bytes(p_size_bytes, byteorder='big')
    return c1_bytes + bytes(encrypted_msg)


  def decrypt_msg(self, msg: bytes) -> bytes:
    # Error handling.
    if self.A_ == None:
      raise Exception("source's factor (A_) not set")
    if self.B_ == None:
      raise Exception("target's factor (B_) not set")

    # Preprocess into byte string of valid size.
    if msg.__class__ != bytes:
      raise Exception("msg must an instance of bytes")
    # Initialize decryption buffer.
    p_size_bytes = math.ceil(math.log2(p) / 8)
    num_encrypted_blocks = len(msg) // self.block_size_
    num_decrypted_bytes = num_encrypted_blocks * self.block_size_
    decrypted_msg = bytearray(num_decrypted_bytes)

    # Decrypt.
    c1_bytes = msg[0:p_size_bytes]
    encrypted_msg = msg[p_size_bytes:]
    c1_int = int.from_bytes(c1_bytes, byteorder='big')
    c1_pow_a_inv = pow(c1_int, self.p_ - 1 - self.private_exp_, self.p_)
    for i in range(0, num_encrypted_blocks):
      encrypted_slice = slice(self.block_size_ * i, self.block_size_ * (i + 1))
      c2_bytes = encrypted_msg[encrypted_slice]
      decrypted_msg[encrypted_slice] = self._decrypt_block(c2_bytes, c1_pow_a_inv)

    return bytes(decrypted_msg)
