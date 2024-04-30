<script setup>
const taskId = ref(null);
const props = defineProps(['url']);
async function runJob() {
    try{
        const response = await axios.get(props.url);
        taskId.value = response.data;
    } catch(error) {
        alert(error);
    }
}
</script>

<template>
    <div>
        <button @click="runJob">run</button>
        <a :href="'/download/' + taskId">{{ taskId }}</a>
    </div>
</template>