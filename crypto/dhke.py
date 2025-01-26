import random
import math
import secrets

class DiffieHellmanClient:
  def __init__(self, **kwargs):
    """
    Client class for executing the Diffie-Hellman key exchange. 
    
    Implements encryption and decryption in the Elgamal cryptosystem. Does not 
    implement any protocols for network communication.
    """
    # Bytes per block.
    self._block_size = 128
    if "block_size" in kwargs:
      self._block_size = kwargs["block_size"]
    # Mod (`p`: public).
    self._p = None
    if "public_prime" in kwargs:
      self._p = kwargs["public_prime"]
    # Base (`g`: public).
    self._g = None
    if "public_base" in kwargs:
      self._g = kwargs["public_base"]
    # Exponent (`a`: private) and my shared factor (`A` = `g`^`a`: public).
    self._private_exp = None
    self._A = None
    if kwargs.get("auto_set_private_exp", False):
      self.set_private_exp()
      self.compute_g_pow_a()
    # Your shared factor (`B` = `g`^`b`: public).
    self._B = None
    # Shared public key ((`g`^`a`)^`b`: public).
    self._shared_public_key = None

  def block_size(self):
    return self._block_size
  
  def set_public_prime(self, p: int):
    if math.ceil(math.log2(p) / 8) < self._block_size:
      raise Exception("prime size (in bytes) must be at least as large as block size.")
    self._p = p

  def set_public_base(self, g: int):
    self._g = g

  def set_private_exp(self):
    if self._p == None:
      raise Exception("mod (p_) not set")
    self._private_exp = secrets.randbelow(self._p - 1) + 1

  def compute_g_pow_a(self):
    if self._p == None:
      raise Exception("mod (p_) not set")
    if self._g == None:
      raise Exception("base (g_) not set")
    if self._private_exp == None:
      self.set_private_exp()
    self._A = pow(self._g, self._private_exp, self._p)

  def set_shared_public_key(self, B: int):
    if self._p == None:
      raise Exception("mod (p_) not set")
    if self._g == None:
      raise Exception("base (g_) not set")
    if self._private_exp == None:
      self.set_private_exp()
    if self._A == None:
      self.compute_shared_public_key_factor()
    self._B = B
    self._shared_public_key = pow(B, self._private_exp, self._p)

  def _decrypt_block(self, msg_block: bytes, c1_pow_a_inv: int) -> bytes:
    c2_int = int.from_bytes(msg_block)
    decrypted_block_int = (c1_pow_a_inv * c2_int) % self._p
    decrypted_block_bytes = decrypted_block_int.to_bytes(self._block_size, byteorder='big')
    return decrypted_block_bytes

  def _encrypt_block(self, encrypted_block: bytes, k: int) -> bytes:
    encrypted_block_int = int.from_bytes(encrypted_block)
    c2_int = (encrypted_block_int * pow(self._B, k, self._p)) % self._p
    c2_bytes = c2_int.to_bytes(self._block_size, byteorder='big')
    return c2_bytes

  def encrypt_msg(self, msg: str | bytes) -> bytes:
    # Error handling.
    if self._A == None:
      raise Exception("source's factor (A_) not set")
    if self._B == None:
      raise Exception("target's factor (B_) not set")

    # Preprocess into byte string of valid size.
    if msg.__class__ not in (str, bytes):
      raise Exception("msg must be str or bytes")
    if msg.__class__ == str:
      msg = bytes(msg, encoding='utf-8') # Always UTF 8
    p_size_bytes = math.ceil(math.log2(p) / 8)
    num_encrypted_blocks = math.ceil(len(msg) / self._block_size)
    num_encrypted_bytes = num_encrypted_blocks * self._block_size
    encrypted_msg = bytearray(num_encrypted_bytes)

    # Generate temporary `k`.
    k = secrets.randbelow(self._p)

    # Encrypt.
    c1_int = pow(self._g, k, self._p)
    for i in range(0, num_encrypted_blocks):
      encrypted_slice = slice(self._block_size * i, self._block_size * (i + 1))
      decrypted_block_bytes = msg[encrypted_slice]
      encrypted_msg[encrypted_slice] = self._encrypt_block(decrypted_block_bytes, k)

    c1_bytes = c1_int.to_bytes(p_size_bytes, byteorder='big')
    return c1_bytes + bytes(encrypted_msg)


  def decrypt_msg(self, msg: bytes) -> bytes:
    # Error handling.
    if self._A == None:
      raise Exception("source's factor (A_) not set")
    if self._B == None:
      raise Exception("target's factor (B_) not set")

    # Preprocess into byte string of valid size.
    if msg.__class__ != bytes:
      raise Exception("msg must an instance of bytes")
    # Initialize decryption buffer.
    p_size_bytes = math.ceil(math.log2(p) / 8)
    num_encrypted_blocks = len(msg) // self._block_size
    num_decrypted_bytes = num_encrypted_blocks * self._block_size
    decrypted_msg = bytearray(num_decrypted_bytes)

    # Decrypt.
    c1_bytes = msg[0:p_size_bytes]
    encrypted_msg = msg[p_size_bytes:]
    c1_int = int.from_bytes(c1_bytes, byteorder='big')
    c1_pow_a_inv = pow(c1_int, self._p - 1 - self._private_exp, self._p)
    for i in range(0, num_encrypted_blocks):
      encrypted_slice = slice(self._block_size * i, self._block_size * (i + 1))
      c2_bytes = encrypted_msg[encrypted_slice]
      decrypted_msg[encrypted_slice] = self._decrypt_block(c2_bytes, c1_pow_a_inv)

    return bytes(decrypted_msg)
