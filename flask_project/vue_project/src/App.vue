<script setup>
import AlignReads from './AlignReads.vue'
</script>

<template>
  <div class="outer">
    <h1>Align Reads</h1>
    <AlignReads />
  </div>

  <div class="outer">
    <h1>Monitor Progresses</h1>
    <a href="/flower" target="_blank">Celery Flower</a>
  </div>
</template>

<style>
    .outer {
       outline: 0px solid #CCC;
       border: 50px solid #FFF;
       background-color: #999;
    }
</style>
